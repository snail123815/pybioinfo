from dataclasses import dataclass


@dataclass
class HmmerHomologousProtConfig:
    MIN_PROTEIN_LEN: int = 60
    T_E: float = (
        1e-5
        # -E, *report* sequences <= E in output
    )
    T_INCE: float = (
        1e-10
        # --incE, *significance*, consider seqs <= incE as significance
    )
    T_DOME: float = (
        1e-5
        # --domE, *report* domains <= domE
    )
    T_INCDOME: float = (
        1e-10
        # --incdomE, *significance*, consider domains <= incdomE as significance
    )
    GATHER_T_E: float = (
        1e-10
        # Domtbl filter, on column 7:
        # E-value of the overall sequence/profile comparison (including all domains)
    )
    GATHER_T_DOME: float = (
        1e-20
        # Used in cal_cov(), coverage calculation
        # Only domE <= this value will be considered
    )
    GATHER_T_COV: float = 0.7  # Coverage of both query and target
    LEN_DIFF: float = 0.2  # diff / min(qlen, tlen)


hhpc = HmmerHomologousProtConfig()
