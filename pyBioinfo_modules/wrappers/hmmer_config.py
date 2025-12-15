from dataclasses import dataclass


@dataclass
class BaseHmmerConfig:
    """
    Default configuration for HMMER tools.
    """

    T_E: float = 10  # -E, *report* sequences <= E in output
    T_INCE: float = 0.001  # --incE, consider seqs <= incE as *significance*
    T_DOME: float = 10  # --domE, *report* domains <= domE
    T_INCDOME: float = (
        0.001  # --incdomE, consider domains <= incdomE as *significance*
    )


@dataclass
class HmmerHomologousProtConfig(BaseHmmerConfig):
    """
    hhpc = Hmmer Homologous Proteins search Config
    Configuration for homologous protein search using HMMER.
    """

    T_E: float = 1e-5
    T_INCE: float = 1e-10
    T_DOME: float = 1e-5
    T_INCDOME: float = 1e-10
    MIN_PROTEIN_LEN: int = 60
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
