import numpy as np
import numpy.polynomial.polynomial as poly
import pandas as pd


def filter_peaks(peak_df, filter_groups: list[list[str, list]]) -> pd.DataFrame:
    working_peak_df = peak_df.copy()
    for filter_group in filter_groups:
        filter_key, filter_limits = filter_group
        if filter_key is None:
            continue
        if filter_key not in peak_df.columns:
            print(
                f"Filter key {filter_key} not found in peak table, ignored "
                f"filtering {filter_group}"
            )
        if filter_key.lower() == "chr":
            # Maybe generalize to str type filters
            for target_chr in filter_limits:
                working_peak_df = working_peak_df.query(
                    f"chr == '{target_chr}'"
                )
        else:
            working_peak_df = working_peak_df.query(
                f"{filter_key} >= {filter_limits[0]} and "
                f"{filter_key} <= {filter_limits[1]}"
            )
    return working_peak_df


def findPeakTrough(array):
    """findPeakTrough(array) -> min_indices, max_indices
    Finds the turning points within an 1D array and returns the indices of the minimum and
    maximum turning points in two separate lists.
    https://stackoverflow.com/a/48360671/6823079
    answered Jan 20 2018 at 20:18
    by Timo
    """
    idx_max, idx_min = [], []
    if len(array) < 3:
        return idx_min, idx_max

    NEUTRAL, RISING, FALLING = range(3)

    def get_state(a, b):
        if a < b:
            return RISING
        if a > b:
            return FALLING
        return NEUTRAL

    ps = get_state(array[0], array[1])
    begin = 1
    for i in range(2, len(array)):
        s = get_state(array[i - 1], array[i])
        if s != NEUTRAL:
            if ps != NEUTRAL and ps != s:
                if s == FALLING:
                    idx_max.append((begin + i - 1) // 2)
                else:
                    idx_min.append((begin + i - 1) // 2)
            begin = i
            ps = s
    return idx_min, idx_max


def peakLengthAtLoc(data, series):
    if not hasattr(series, "__iter__"):
        series = [
            series,
        ]
    lengthNums = []
    for x in series:
        if 0 <= x <= len(data.index):
            lengthNums.append(data.length.iloc[np.int(x)])
        else:
            lengthNums.append("")
    if len(lengthNums) == 1:
        return lengthNums[0]
    return lengthNums


def minMaxPolyfit(x, y, degree=10):
    """Fit a polynomial curve to data and get
    the turning points (peak and trough) on this curve.
    x,y,degree = 10
    x,y are list of values"""
    coefs = poly.polyfit(x, y, degree)
    Fit = poly.Polynomial(coefs)
    yTrend = [
        np.NaN,
    ]
    yTrend_trend = [
        np.NaN,
    ]
    for index in range(len(y) - 1):
        if index < 500:
            yTrend.append(np.NaN)
            yTrend_trend.append(np.NaN)
        else:
            yTrend.append(Fit(index + 1) - Fit(index))
            if np.isnan(yTrend[index - 1]):
                yTrend_trend.append(np.NaN)
            else:
                # this is to find real turning point - tangent angle change
                yTrend_trend.append(yTrend[index] - yTrend[index - 1])
    min_indices, max_indices = findPeakTrough(yTrend_trend)
    return min_indices, max_indices, Fit


def filterFoldEnrichment(peakDF, method=["single", 10]):
    methodIsList = False
    if isinstance(method, list):
        methodIsList = True
        thresh = method[1]
        method = method[0]
    else:
        thresh = 10

    if method not in ["single", "max", "min"]:
        raise Exception(f"Wrong method '{method}' ('single','max','min')")
    if method == "single":
        if "fold_enrichment" not in peakDF.columns:
            raise Exception(
                f"Wrong dataframe for {method} fold enrichment filter: {peakDF.columns}"
            )
        y = peakDF.fold_enrichment
        if len(y) <= 2000 or methodIsList:
            print(f"No polyfit, setting thresh for fold_enrichment = {thresh}")
        else:
            bondarieMin = 1300
            x = range(1, len(y) + 1)
            min_indices, max_indices, Fit = minMaxPolyfit(x, y)
            thresh = min(
                [
                    Fit(lim)
                    for lim in [
                        indice
                        for indice in max_indices
                        if indice >= bondarieMin
                    ]
                ]
            )
            print(f"Done polyfit, thresh for fold_enrichment = {thresh}")
    if method == "max":
        if "fold_enrichment_A" not in peakDF.columns:
            raise Exception(
                f"Wrong dataframe for {method} fold enrichment filter: {peakDF.columns}"
            )
        maxFold = peakDF.loc[:, ["fold_enrichment_A", "fold_enrichment_B"]].max(
            axis=1
        )
        peakDF = peakDF.assign(fold_enrichment=maxFold)
    if method == "min":
        if "fold_enrichment_A" not in peakDF.columns:
            raise Exception(
                f"Wrong dataframe for {method} fold enrichment filter: {peakDF.columns}"
            )
        minFold = peakDF.loc[:, ["fold_enrichment_A", "fold_enrichment_B"]].min(
            axis=1
        )
        peakDF = peakDF.assign(fold_enrichment=minFold)
    filtered = peakDF.query("fold_enrichment >= @thresh")
    return filtered


# filterFoldEnrichment method = ['single','max','min']


def filterLength(peakDF, method="dist"):
    peakDF = peakDF.assign(length=peakDF.end - peakDF.start)
    peakDF.sort_values("length", inplace=True)
    if method == "dist":
        medians = peakDF.length.median()
        minLength = int(peakDF.length.quantile(q=[0.25]))
        maxLength = int(peakDF.length.quantile(q=[0.75]))
    elif method == "polyfit":
        # data has to be origional peak files
        # this has to be checked by plotting in peak_enrichment_modeling
        y = peakDF.fold_enrichment
        if len(y) <= 2000:
            raise Exception(
                "You are doing polyfit in a data set less then 2000 entries."
            )
        bondarieMin = 1300
        x = range(1, len(y) + 1)
        min_indices, max_indices, Fit = minMaxPolyfit(x, y)
        foldThresh = min(
            [
                Fit(lim)
                for lim in [
                    indice for indice in max_indices if indice >= bondarieMin
                ]
            ]
        )
        max_indices = [
            indice for indice in max_indices if indice >= bondarieMin
        ]
        if len(max_indices) != 2:
            raise Exception(
                f"Please check the ployfit graph, there are {len(max_indices)} turning points for now."
            )
        # gether data for output
        minLength, maxLength = [
            peakLengthAtLoc(peakDF, indice)
            for indice in max_indices
            if indice >= bondarieMin
        ]
    elif isinstance(method, list):
        try:
            minLength, maxLength = method
        except BaseException:
            raise Exception(
                "Method should be [minLength, maxLength] if you do not fillter"
            )
    else:
        raise Exception(f"Method {method} wrong for length filter")
    filtered = peakDF.query("@minLength <= length <= @maxLength")
    filtered.sort_index()
    return filtered


# filterLength method = ['dist','polyfit']


def change_location_to_summit(
    peak_df, to_left=150, to_right=150, seq_end_at=None
):
    if "abs_summit" not in peak_df.columns:
        raise Exception("No summit info in dataframe")
    peak_range_max = peak_df.end.max()
    peak_df.start = peak_df.abs_summit - to_left
    peak_df.end = peak_df.abs_summit + to_right
    if seq_end_at:
        peak_df.end = min(peak_df.end, seq_end_at)
    else:
        peak_df.end = peak_df.end.clip(upper=peak_range_max)
    peak_df.start = peak_df.start.clip(lower=0)
    return peak_df


def filterLikely(peakDF, method):
    thresh = float(method)
    if "likely_difference" in peakDF.columns:
        if thresh >= 10:
            raise Exception(
                "Threshold error, for common peaks the value should be the lesser the better, set up a small value ~1"
            )
        filtered = peakDF.query("likely_difference <= @thresh")
    elif "log10_likely" in peakDF.columns:
        if thresh < 10:
            raise Exception(
                "Threshold error, for unique peaks, the value should be the bigger the better, set up a big value >=10"
            )
        filtered = peakDF.query("log10_likely >= @thresh")
    elif "-log10(pvalue)" in peakDF.columns:
        if thresh < 10:
            raise Exception(
                "Threshold error, for peaks, the value should be the bigger the better, set up a big value >=10"
            )
        filtered = peakDF.query("`-log10(pvalue)` >= @thresh")

    else:
        raise Exception(
            f"Likely filter failed for this data frame with header\n{peakDF.columns}"
        )
    return filtered


# filterLikely
