import subprocess
from os.path import commonpath
from pathlib import Path
from tempfile import NamedTemporaryFile, _TemporaryFileWrapper

IMPLEMENTED_COMPRESSION_FORMATS: list[str] = [".gz", ".xz"]


def compressFile(
    filePath: Path, compressionFormat: str, keepOrigion=False, dry=False
) -> Path:
    match compressionFormat:
        case ".gz":
            prog = "gzip"
        case ".xz":
            prog = "xz"
        case _:
            raise NotImplementedError(
                f"Decompress {filePath.suffix} file is not supported."
            )
    assert filePath.suffix != compressionFormat
    resultFile = filePath.parent / (filePath.name + compressionFormat)
    if dry:
        return resultFile

    with open(resultFile, "wb") as rf:
        compress = subprocess.run(
            [prog, "-c", filePath.resolve()], stdout=rf, stderr=subprocess.PIPE
        )
        assert compress.returncode == 0

    if not keepOrigion:
        assert resultFile.exists()
        filePath.unlink()
    return resultFile


def splitStemSuffixIfCompressed(
    filePath: Path,
    allowedFormats: list[str] = IMPLEMENTED_COMPRESSION_FORMATS,
    fullSuffix: bool = False,
) -> tuple[str, str]:
    if filePath.suffix in allowedFormats:
        stem = filePath.with_suffix("").stem
        suffix = filePath.with_suffix("").suffix
        if fullSuffix:
            suffix = filePath.suffix + suffix
        return stem, suffix
    else:
        return filePath.stem, filePath.suffix


def getSuffixIfCompressed(
    filePath: Path,
    allowedFormats: list[str] = IMPLEMENTED_COMPRESSION_FORMATS,
    fullSuffix: bool = False,
) -> str:
    return splitStemSuffixIfCompressed(filePath, allowedFormats, fullSuffix)[1]


def getStemIfCompressed(
    filePath: Path, allowedFormats: list[str] = IMPLEMENTED_COMPRESSION_FORMATS
) -> str:
    return splitStemSuffixIfCompressed(filePath, allowedFormats)[0]


def decompressFile(filePath: Path) -> Path:
    match filePath.suffix:
        case ".gz":
            prog = "gzip"
        case ".xz":
            prog = "xz"
        case _:
            raise NotImplementedError(
                f"Decompress {filePath.suffix} file is not supported."
            )
    resultFilePath = filePath.with_suffix("")
    with open(resultFilePath, "w") as rf:
        decompress = subprocess.run(
            [prog, "-dc", filePath.resolve()], stdout=rf, stderr=subprocess.PIPE
        )
    if not resultFilePath.exists():
        raise FileNotFoundError(
            "\n".join(
                [
                    f"Unzip file {filePath} failed: no output file found:",
                    " ".join(decompress.args),
                    decompress.stderr.decode(),
                ]
            )
        )
    return resultFilePath


def decompFileIfCompressed(
    filePath: Path, allowedFormats: list[str] = IMPLEMENTED_COMPRESSION_FORMATS
) -> tuple[Path, bool]:
    if filePath.suffix in allowedFormats:
        return decompressFile(filePath), True
    else:
        return filePath, False


def decompressToTempTxt(filePath: Path) -> _TemporaryFileWrapper:
    # To access the file outside this function,
    # do NOT: a = decompressToTempTxt(filePath).name
    # it will kill the temporary file.
    # You can do with decompressToTempTxt(filePath) as t:
    match filePath.suffix:
        case ".gz":
            prog = "gzip"
        case ".xz":
            prog = "xz"
        case _:
            raise NotImplementedError(
                f"Decompress {filePath.suffix} file is not supported."
            )
    outTempFile = NamedTemporaryFile()
    with open(outTempFile.name, "w") as out:
        subprocess.run([prog, "-dc", filePath], stdout=out, check=True)
    return outTempFile


def getRootAndFiles(
    pathList: list[Path], allowedExts
) -> tuple[Path, list[Path]]:
    targetFiles: list[Path] = []
    for pathIn in pathList:
        pathIn = Path(pathIn).resolve()
        if pathIn.is_dir():
            for f in pathIn.iterdir():
                if getSuffixIfCompressed(f) in allowedExts:
                    targetFiles.append(f)
        elif pathIn.is_file():
            assert getSuffixIfCompressed(pathIn) in allowedExts
            targetFiles.append(pathIn)
    rootPath: Path = Path(commonpath(targetFiles))
    if rootPath.is_file():
        rootPath = rootPath.parent
    return rootPath, targetFiles
