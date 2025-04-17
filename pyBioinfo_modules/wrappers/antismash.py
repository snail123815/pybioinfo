import logging
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path, PurePath
from typing import Literal, TypedDict

import ijson
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from pyBioinfo_modules.basic.decompress import decompFileIfCompressed
from pyBioinfo_modules.bio_sequences.bio_seq_file_extensions import \
    FNA_EXTENSIONS
from pyBioinfo_modules.wrappers._environment_settings import (
    ANTISMASH_ENV, CONDAEXE, SHELL, withActivateEnvCmd)

logger = logging.getLogger(__name__)


def findClusterNumberStr(file: Path, numberOnly: bool = False) -> str | None:
    match = clusterNumberPattern.search(file.name)
    if match:
        if not numberOnly:
            return match[0].split(".")[1]
        else:
            return match[0].split(".")[1][-3:]
    else:
        return None


antismashClusterGbkFileNameTest = "NNNNNNNNNNNn.region001.gbk"
clusterGbkGlobTxt = r"*region[0-9][0-9][0-9].gbk"
assert PurePath(antismashClusterGbkFileNameTest).match(clusterGbkGlobTxt)
clusterNumberPattern = re.compile(r"\.region[0-9]{3}\.gbk$")
assert (
    findClusterNumberStr(Path(antismashClusterGbkFileNameTest)) == "region001"
)


def log_antismash_version(
    condaexe=CONDAEXE, antismash_env=ANTISMASH_ENV, shell=SHELL
) -> str:
    logger.info(f"antiSMASH environment location: {antismash_env}")
    logger.info(f"Conda executable: { condaexe }")
    logger.info(f"Shell: { shell }")
    cmd = withActivateEnvCmd(
        "antismash --version", antismash_env, condaexe, shell
    )
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            shell=True,
            executable=shell,
        )
        result.check_returncode()
        version = result.stdout.decode()
    except subprocess.CalledProcessError as e:
        logger.error(f"Command '{cmd}' failed")
        raise e
    version = version.split(" ")[1].strip()
    logger.info(f"antiSMASH version: {version}")
    return version


def runAntismash(
    inputFilePath: Path,
    title: str | None = None,
    description: str | None = None,
    taxon: Literal["bacteria", "fungi"] = "bacteria",
    completeness: Literal[1, 2, 10] = 2,
    condaExe: Literal["conda", "mamba", "micromamba"] = CONDAEXE,
    condaEnv: Path | None = ANTISMASH_ENV,
    cpu: int = 4,
    output: Path | None = None,
    shell: Literal["bash", "zsh"] = SHELL,
    outputPrefix: str = "antismash",
    addDateTimeToPrefix: bool = False,
    geneFinding: Literal[
        "glimmerhmm", "prodigal", "prodigal-m", "none", "error"
    ] = "none",
    # There will be short contig without CDS, if set to 'error'
    # program will directly stop.
    defaultGeneFinding: str = "prodigal",
    silent: bool = False,
    dry: bool = False,
    overwrite: bool = False,
    existsOk: bool = False,
) -> Path:

    logger.info(f"Running antiSMASH for {inputFilePath}")
    log_antismash_version(condaExe, condaEnv, shell)

    inputFilePath, unzip = decompFileIfCompressed(inputFilePath)

    outputPrefix = "_".join(
        item
        for item in [
            outputPrefix,
            title,
            f"level{completeness}",
        ]
        if item is not None
    )
    if addDateTimeToPrefix:
        timeStr = datetime.now().strftime(r"%Y%m%d%H%M")
        outputPrefix += "_" + timeStr
    try:
        if output is None:
            outdir = inputFilePath.parent / outputPrefix
        elif (
            output.is_dir()
            and output.resolve() in inputFilePath.resolve().parents
        ):
            # If output is a parent of inputfile, also create an output dir.
            outdir = output / outputPrefix
        else:
            outdir = output

        if (outdir / "index.html").exists():
            if overwrite:
                shutil.rmtree(outdir)
            elif existsOk:
                logger.info(f"Find result file in {outdir}, pass.")
                return outdir.resolve()
            else:
                raise FileExistsError(str(outdir))
        elif outdir.exists():
            shutil.rmtree(outdir)

        cmd = (
            f"antismash --cpus {cpu}"
            + " --minimal"
            + " --skip-zip-file"
            + f" --taxon {taxon}"
            + f" --html-title {outputPrefix}"
            + f" --output-dir {outdir}"
        )
        if inputFilePath.suffix in FNA_EXTENSIONS and geneFinding == "auto":
            cmd += f" --genefinding-tool {defaultGeneFinding}"
        elif geneFinding == "auto":
            cmd += f" --genefinding-tool none"
        else:
            cmd += f" --genefinding-tool {geneFinding}"
        if description is not None:
            cmd += f" --html-description {description}"

        if completeness >= 2:
            cmd = cmd.replace(" --minimal", "")
            cmd += " --cb-knownclusters"
            cmd += " --cb-subclusters"
            cmd += " --asf"
        if completeness >= 3:
            cmd += " --cb-general"
            cmd += " --cc-mibig"
            cmd += " --clusterhmmer"
            cmd += " --pfam2go"
            if taxon == "fungi":
                cmd += " --cassis"
        if completeness >= 4:
            cmd += " --rre"
            cmd += " --fullhmmer"
            cmd += " --tigrfam"
            cmd += " --smcog-trees"

        cmd += f" {inputFilePath}"

        if not silent:
            logger.info(cmd)

        cmd = withActivateEnvCmd(cmd, condaEnv, condaExe, shell)

        if dry:
            logger.info(cmd)
        else:
            commandResult = subprocess.run(
                cmd, capture_output=True, shell=True, executable=shell
            )
            if commandResult.returncode != 0:
                logger.error("Failed antismash:")
                logger.error(cmd)
                logger.error(commandResult.stdout.decode())
                logger.error(commandResult.stderr.decode())
    finally:
        if unzip:
            os.remove(str(inputFilePath))
    logger.info(f"Done antiSMASH for {inputFilePath}")

    return outdir.resolve()


class ClusterInfo(TypedDict):
    gbkFile: Path
    gcProducts: str
    organism: str
    fromSequence: str
    coreRelativeLocs: list[FeatureLocation]
    joinedProteinFastaFile: Path
    fastaId: str
    # TODO
    # Each AntiSMASH cluster is a "region", feature type == 'region',
    # parse region info in 'qualifiers'
    #     "region_number"
    #     "product" (could be multiple)
    # Each AntiSMASH cluster might have multiple protocluster,
    # feature type == 'protocluster',
    # parse protocluster info in 'qualifiers'
    #     "protocluster_number"
    #     "product"


def parseClusterGbk(
    infile: Path,
    proteinFastaDir: Path,
    id: int | None = None,
    nflank: int = 0,
) -> ClusterInfo:
    """Parses the genbank files for DNA, protein, cluster, organism
    [Rewrote of parsegbkcluster() from BiG-MAP]
    parameters
    ----------
    infile
        Path of .gbk file
    proteinFastaDir
        Path of protein fastas to store joined proteins for each cluster
    nflank
        Number of CDS to flank the core regions
    returns
    ----------
    DNA = DNA sequence
    proteins = protein sequence
    clustername = name of the cluster
    organism = name of the organism
    #cluster_enzymes = {loc:gene_kind}
    """
    proteins = []
    gcProducts = []
    cdsIndexs = []
    coreIndexs = []
    coreRelativeLocs = []
    records = list(SeqIO.parse(infile, "genbank"))
    assert len(records) == 1
    record = records[0]
    for i, feature in enumerate(record.features):
        # Parsing the regionname, part of antismash output
        if "region" in feature.type:
            if "product" in feature.qualifiers:
                for product in feature.qualifiers["product"]:
                    gcProducts.append(product)
        # Parsing the protein sequence
        if feature.type == "CDS":
            cdsIndexs.append(i)
            # remembering every CDS gene index
            try:
                proteins.append(feature.qualifiers["translation"][0])
            except KeyError:
                pass
            # Parsing the relative core locations
            try:
                if feature.qualifiers["gene_kind"][0] == "biosynthetic":
                    coreIndexs.append(i)
            except KeyError:
                pass

    if len(coreIndexs) > 0:
        if cdsIndexs.index(min(coreIndexs)) - nflank < 0 or cdsIndexs.index(
            max(coreIndexs)
        ) + nflank + 1 > len(cdsIndexs):
            print(
                "!!!flank_genes (-f) is higher than the number of flanking "
                + f"genes in the cluster of file: {infile}, "
                + "using whole gene cluster instead!!!"
            )
        if nflank == 0:
            coreRelativeLocs = [record.features[i].location for i in coreIndexs]
        else:
            if cdsIndexs.index(min(coreIndexs)) - nflank >= 0:
                core_region = cdsIndexs[
                    cdsIndexs.index(min(coreIndexs))
                    - nflank : cdsIndexs.index(max(coreIndexs))
                    + nflank
                    + 1
                ]
            else:
                core_region = cdsIndexs[
                    0 : cdsIndexs.index(max(coreIndexs)) + nflank + 1
                ]
            coreRelativeLocs = [
                record.features[i].location for i in core_region
            ]

    organism = record.description
    # Parsing organism name
    # The attr. description is from "DEFINATION" line for gbk record
    # biopython/Bio/GenBank/__init__.py line 764
    # in class _FeatureConsumer(_BaseGenBankConsumer)
    if " " in organism:
        organism = "_".join(organism.split(",")[0].split())
    organism = re.sub(f"strain", "", organism)
    organism = re.sub(f'[()"]+', "", organism)
    organism = re.sub(f"[\\/]+", "-", organism)
    organism = re.sub(f"[.]+", ".", organism)
    organism = re.sub(f"[_]+", "_", organism)

    # write to protein fasta file
    joinProteins = Seq("".join(proteins))
    regionAndNumber = findClusterNumberStr(infile)
    fromSequence = infile.name.split(regionAndNumber)[0][:-1]
    if id is None:
        proteinFastaFileName = (
            f"GC_PROT-{fromSequence}-{regionAndNumber}" + ".fasta"
        )
    else:
        proteinFastaFileName = f"P{id}"
    fastaId = (
        f"{organism}|{fromSequence}|"
        + f"GC_PROT--{regionAndNumber}"
        + f"--Entryname={':'.join(gcProducts)}"
    )
    proteinFastaFile = proteinFastaDir / proteinFastaFileName
    SeqIO.write(
        SeqRecord(joinProteins, id=fastaId, description=""),
        proteinFastaFile,
        "fasta",
    )
    assert proteinFastaFile.is_file()

    return ClusterInfo(
        gbkFile=infile,
        gcProducts=":".join(gcProducts),
        organism=organism,
        fromSequence=fromSequence,
        coreRelativeLocs=coreRelativeLocs,
        joinedProteinFastaFile=proteinFastaFile,
        fastaId=fastaId,
    )


def find_NRPS_TE_domain(jsonResultPath: Path) -> list[SeqRecord]:
    teDomains: list[SeqRecord] = []
    resultName = jsonResultPath.stem

    def TE_follow_PCP(domainList: list[str]) -> bool:
        try:
            locationTe = domainList.index("Thioesterase")
            try:
                if domainList[locationTe - 1] == "PCP":
                    return True
            except IndexError:
                pass
            try:
                if domainList[locationTe + 1] == "PCP":
                    return True
            except IndexError:
                pass
        except ValueError:
            pass
        return False

    with open(jsonResultPath, "rb") as jsonHandle:
        asResultRecords = ijson.kvitems(jsonHandle, "records.item")
        recordFeatures = (v for k, v in asResultRecords if k == "features")
        recordIds = (v for k, v in asResultRecords if k == "id")

        for id, rec in zip(recordIds, recordFeatures):
            regionNumber = -1
            protoclusterNumber = -1
            protoclusterProduct = ""
            domainList: list[str] = []

            for feat in rec:
                if feat["type"] == "region":
                    regionNumber = int(feat["qualifiers"]["region_number"][0])
                elif feat["type"] == "protocluster":
                    protoclusterNumber = int(
                        feat["qualifiers"]["protocluster_number"][0]
                    )
                    protoclusterProduct = feat["qualifiers"]["product"][0]
                if (
                    protoclusterNumber == -1
                    or not protoclusterProduct == "NRPS"
                ):
                    continue

                if feat["type"] == "aSDomain":
                    domain = feat["qualifiers"]["aSDomain"][0]
                    domainList.append(domain)
                    if domain == "Thioesterase":
                        teSeqRec = SeqRecord(
                            Seq(feat["qualifiers"]["translation"][0]),
                            id="_".join(
                                feat["qualifiers"]["domain_id"][0].split("_")[
                                    1:
                                ]
                            ),
                            name=feat["qualifiers"]["locus_tag"][0],
                            description=(
                                f"    {resultName}_seq_{id}"
                                f"_region{regionNumber:0>3}"
                                f"_protocluster{protoclusterNumber:0>3}"
                            ),
                        )

                # See if next CDS started
                if feat["type"] == "CDS":
                    if TE_follow_PCP(domainList):
                        assert "teSeqRec" in locals()
                        tqdm.write(
                            "Found valid TE domain:\n"
                            f"    {resultName}_seq_{id}"
                            f"_region{regionNumber:0>3}"
                            f"_protocluster{protoclusterNumber:0>3}"
                        )
                        teDomains.append(teSeqRec)
                    domainList = []

    jsonHandle.close()
    return teDomains
