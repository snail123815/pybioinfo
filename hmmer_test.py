from pathlib import Path

from pyBioinfo_modules.wrappers.hmmer import (
    full_proteome_jackhmmer,
    run_jackhmmer,
)
from pyBioinfo_modules.wrappers.hmmer_config import (
    BaseHmmerConfig,
    HmmerHomologousProtConfig,
)

with open(Path("/home/duc/more.fa"), "r") as f:
    a = run_jackhmmer(
        query_fasta=f,
        target_fasta_path=Path(
            "/home/duc/uniprotkb_streptomyces_coelicolor_a3_AN_2025_04_22.fasta"
        ),
        domtblout_path=Path("/home/duc/more.direct_read.domtblout"),
        jackhmmer_cfg=BaseHmmerConfig(T_DOME=1e-10),
    )

print(a)

full_proteome_jackhmmer(
    query_fasta_path=Path("/home/duc/more.fa"),
    target_fasta_path=Path(
        "/home/duc/uniprotkb_streptomyces_coelicolor_a3_AN_2025_04_22.fasta"
    ),
    domtblout_path=Path("/home/duc/more.full.domtblout"),
    prots_per_group=2,
    hhpc=HmmerHomologousProtConfig(),
)

with open("/home/duc/more.full.domtblout", "r") as f:
    for line in f:
        print(line.strip())
