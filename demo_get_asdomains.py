from pyBioinfo_modules.wrappers.antismash import get_asdomains_json
from pathlib import Path


targetJson = Path("GCF_008634025.1_ASM863402v1_genomic.json")

get_asdomains_json(targetJson)