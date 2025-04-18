from dataclasses import dataclass


@dataclass
class HmmerHomologousProtConfig:
    MIN_PROTEIN_LEN: int = 60
    T_E: float = 1e-5
    T_INCE: float = 1e-10
    T_DOME: float = 1e-5
    T_INCDOME: float = 1e-10
    GATHER_T_E: float = 1e-10
    GATHER_T_DOME: float = 1e-20
    GATHER_T_COV: float = 0.7
    LEN_DIFF: float = 0.2  # diff / min(qlen, tlen)


hhpc = HmmerHomologousProtConfig()
