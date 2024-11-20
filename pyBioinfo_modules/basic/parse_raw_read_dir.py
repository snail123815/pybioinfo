import logging
from pathlib import Path

from pyBioinfo_modules.basic.decompress import splitStemSuffixIfCompressed

logger = logging.getLogger(__name__)


def imputePeSuffix(rawFileNames: list[str], peSfx: list[str] = []):
    if len(peSfx) != 0:
        return peSfx
    assert len(rawFileNames) % 2 == 0 and len(rawFileNames) != 0, (
        "Pair end reads should be in pairs, however"
        f" {len(rawFileNames)} found."
    )
    peFound = False
    for i in range(1, 12):
        s: set[str] = set(fn[-i:] for fn in rawFileNames)
        if len(s) == 2:
            peSfx = sorted(s)
            peFound = True
        if len(s) != 2 and peFound:
            break
    if peFound:
        logger.info(f"Found pairend suffix {peSfx}")
        return peSfx
    raise ValueError(f"pair end suffix not found in {rawFileNames}")


def get_read_files_per_sample(
    parent_dir: list[Path],
    sample_names: list[str] | None = None,
    isPe: bool = False,
) -> tuple[dict[str, list[Path]], list[str], str]:
    # Gether all files, max depth 2
    file_paths: list[Path] = []
    for d in parent_dir:
        for n in d.iterdir():
            if n.is_file():
                file_paths.append(n)
            elif n.is_dir():
                file_paths.extend(f for f in n.iterdir() if f.is_file())
    for f in file_paths.copy():
        ext = splitStemSuffixIfCompressed(f)[1]
        if len(ext.split(".")) < 2 or ext.split(".")[1] not in ["fastq", "fq"]:
            file_paths.remove(f)
    assert len(file_paths) > 0, f"Files not found in {parent_dir}"
    file_paths = sorted(file_paths)
    file_parts = [
        splitStemSuffixIfCompressed(f, fullSuffix=True) for f in file_paths
    ]
    file_names = [fp[0] for fp in file_parts]
    file_fullexts = [fp[1] for fp in file_parts]
    assert len(set(file_fullexts)) == 1, (
        "All files should have same format, "
        f"yet multiple found {set(file_fullexts)}."
    )
    assert len(set(file_names)) == len(file_names), (
        "File names has duplicates"
        ", maybe from different dirs?\n"
        f"{set([fn for fn in file_names if file_names.count(fn) > 1])}"
    )

    if sample_names is not None:
        samples = sample_names
    else:
        samples = []

    sample_file_dict: dict[str, list[Path]] = {}
    peSfx: list[str] = []
    if isPe:
        peSfx = imputePeSuffix(file_names)
        if samples == []:
            samples = (
                sorted(list(set(f[: -len(peSfx[0])] for f in file_names)))
                if len(samples) == 0
                else samples
            )
            assert (len(file_paths) / len(samples)) % 2 == 0
        for s in samples:
            fns = [file_paths[i] for i, fn in enumerate(file_names) if s in fn]
            assert len(fns) > 0, f"Did not find files for sample {s}"
            assert len(fns) % 2 == 0, (
                f"Number of files for sample {s} is not paired.\n" f"{fns}"
            )
            sample_file_dict[s] = fns
    else:  # single end reads
        if samples == []:
            for fn, fp in zip(file_names, file_paths):
                sample_file_dict[fn] = [fp]
            samples = file_names
        elif len(samples) == 1:
            sample_file_dict[samples[0]] = file_paths
        else:
            for s in samples:
                fps = [
                    file_paths[i] for i, fn in enumerate(file_names) if s in fn
                ]
                assert len(fps) > 0, f"Did not find files for sample {s}"
                sample_file_dict[s] = fps

    return sample_file_dict, peSfx, file_fullexts[0]
