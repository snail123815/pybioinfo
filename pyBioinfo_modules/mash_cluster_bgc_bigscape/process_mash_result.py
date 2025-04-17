import sys
from pathlib import Path

import numpy as np
from tqdm import tqdm


class mashBGC_ClusteringResult:
    """
    Store major results after mash clustering and medoid calculating.
    Including:
        clusterInfoDict: dict = pickle.load(fh)
        family_distance_matrice: dict[str, list[list[float]]] = pickle.load(fh)
        dict_medoids: dict[str, list[str]] = pickle.load(fh)



    """

    def __init__(
        self,
    ) -> None:
        pass


def calculate_medoid(
    inputDistanceTablePath: Path,  # output file (return) of mashDistance()
    cutOff: float,  # default 0.8
    med: dict[str, list[str]] = {},  # looks like you can pass previous result?
) -> tuple[dict[str, list[str]], dict[str, list[list[float]]]]:
    """
    re-write of function in BiGMAP https://github.com/medema-group/BiG-MAP
    calculates the GCFs based on similarity threshold
    parameters and calculates the medoid of that GCF
    ----------
    outdir
        string, the path to output directory
    cut_off
        float, between 0 and 1
    returns
    ----------
    dict_medoids = {fasta file of medoid: similar fasta files}
    """
    # Parse the input into a dictionary of gene families
    family: dict[str, str] = {}
    # family dict, key = family members, values = family names
    familyFiltered = {}
    family_members: dict[str, list[str]] = {}
    # family_members dict, key = family name, values = list of members
    family_distance_matrices: dict[str, list[list[float]]] = {}
    dict_medoids: dict[str, list[str]] = med

    def add_to_distance_matrix(familyName, refId, queryId, distance):
        def add_new_gene(familyName, idList, id) -> int:
            if id not in idList:
                idList.append(id)
                # Extend distance matrix: One new row, one new column
                for row in family_distance_matrices[familyName]:
                    row.append(0)
                family_distance_matrices[familyName].append([0] * len(idList))
            return idList.index(id)

        index1 = add_new_gene(familyName, family_members[familyName], refId)
        index2 = add_new_gene(familyName, family_members[familyName], queryId)
        family_distance_matrices[familyName][index1][index2] = distance
        family_distance_matrices[familyName][index2][index1] = distance
        return ()

    with inputDistanceTablePath.open("r") as input:
        pbar = tqdm(
            total=inputDistanceTablePath.stat().st_size,
            bar_format=r"{l_bar}{bar}| {n:,.0f}/{total:,.0f} {unit} "
            + r"[{elapsed}<{remaining}, {rate_fmt}{postfix}]",
            unit_scale=1 / 1048576,
            unit="MB",
            desc="Generating families from distance file",
        )
        readSize = 0
        for idx, line in enumerate(input):
            readSize += sys.getsizeof(line) - 50  # length of '\n'
            if idx % 10000 == 0:
                pbar.update(readSize)
                readSize = 0
            if line.startswith("#") or line.strip() == "":
                continue

            # Split into tab-separated elements
            refId, queryId, distanceStr, pValueStr, nHashesStr = (
                line.strip().split("\t")
            )
            sharedNhashes, totalNhashes = (
                int(n) for n in nHashesStr.split("/")
            )
            shareRatio = sharedNhashes / totalNhashes
            distance = float(distanceStr)
            pValue = float(pValueStr)
            # Look up the family of the first gene
            # Each family is named by the first gene of the family
            if queryId in family.keys():
                familyName = family[queryId]
            else:  # init a new family
                familyName = queryId
                family[queryId] = familyName
                family_members[familyName] = []
                family_distance_matrices[familyName] = []
            # filter genes that are already part of another family
            if queryId == family[queryId]:
                familyFiltered[queryId] = family[queryId]
            # If refId and queryId don't overlap, then there are two options:
            # 1. refId and queryId do belong to the same family, and the
            #    "not overlap" is "odd" in this case we accept queryId into
            #    our family
            # 2. refId and queryId actually belong to different families,
            #    and we might have to create that family
            if shareRatio <= cutOff:
                # refId doesn't overlap at all, so put that one into a separate
                # family
                if refId in family.keys():
                    if family[refId] == familyName:
                        # refId is in our family, so record the distance
                        add_to_distance_matrix(
                            familyName, queryId, refId, distance
                        )
                    else:
                        # gene is above cut off or doesn't belong to the family
                        pass
                else:
                    # refId doesn't have a family yet, make one
                    gene1_family_name = refId
                    family[refId] = gene1_family_name
                    family_members[gene1_family_name] = []
                    family_distance_matrices[gene1_family_name] = []
                    # insert refId into that family as only member, with a
                    # distance of 0
                    add_to_distance_matrix(
                        gene1_family_name, refId, refId, distance
                    )
            else:
                # There is some overlap, and we want refId also in this family
                family[refId] = familyName
                add_to_distance_matrix(familyName, queryId, refId, distance)

        pbar.update(readSize)
        pbar.close()
    # For each family: Build a distance matrix, and then work out the medoid
    for familyName in familyFiltered.keys():
        # Calculate the medoid from the distances
        np_array = np.asarray(family_distance_matrices[familyName])
        medoid_index = np.argmin(np_array.sum(axis=0))
        # Create a dictionary using the medoid as key and
        # the family_members as values
        medoidName = family_members[familyName][medoid_index]
        dict_medoids[medoidName] = family_members[familyName]
    return dict_medoids, family_distance_matrices
