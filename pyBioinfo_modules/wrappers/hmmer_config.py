from dataclasses import dataclass


@dataclass
class BaseHmmerTblFilters:
    """
    Default configuration for HMMER tools.
    """

    T_E: float = 10  # -E, *report* sequences <= E in output
    T_INCE: float = 0.001  # --incE, consider seqs <= incE as *significance*
    T_DOME: float = 10  # --domE, *report* domains <= domE
    T_INCDOME: float = (
        0.001  # --incdomE, consider domains <= incdomE as *significance*
    )
    MIN_PROTEIN_LEN: int = 0
    GATHER_T_E: float = (
        10
        # Domtbl filter, on column 7:
        # E-value of the overall sequence/profile comparison (including all domains)
    )
    GATHER_T_DOME: float = (
        1
        # Used in cal_cov(), coverage calculation
        # Only domE <= this value will be considered
    )
    GATHER_T_COV: float = 0.0  # Coverage of both query and target
    LEN_DIFF: float = float("inf")  # diff / min(qlen, tlen)


HmmerHomologousProtFillters = BaseHmmerTblFilters(
    T_E=1e-5,
    T_INCE=1e-10,
    T_DOME=1e-5,
    T_INCDOME=1e-10,
    MIN_PROTEIN_LEN=60,
    GATHER_T_E=1e-10,
    # Domtbl filter, on column 7:
    # E-value of the overall sequence/profile comparison (including all domains)
    GATHER_T_DOME=1e-20,
    # Used in cal_cov(), coverage calculation
    # Only domE <= this value will be considered
    GATHER_T_COV=0.7,
    LEN_DIFF=0.2,
)
