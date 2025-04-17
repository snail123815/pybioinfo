import pickle
from glob import glob
from pathlib import Path

from tqdm import tqdm


def findStartLine(csvFile: Path) -> int:
    with csvFile.open("r") as file:
        for i, line in enumerate(file):
            if not line.strip()[0] == "#":
                if line != "":
                    return i
        return 0


def globFilesSafely(
    sourceDir: Path,
    globPattern: str,
    resultPickle: Path | None = None,
    reload: bool = False,
    showProgress: bool = False,
) -> list[Path]:
    """Make sure glob files will follow symlinks
    In pathlib.Path, the glob() function do not follow symlink:
    [pathlib.Path.glob does not follow symlinks · Issue #77609 · python/cpython
    (github.com)](https://github.com/python/cpython/issues/77609)
    But this is hardly a bug because this is a system behaviour.
    However this will happily follow symlinks:
    ```python
    from glob import glob
    glob(str, recursive=True)
    ```
    """
    if resultPickle is not None and resultPickle.is_file() and not reload:
        with resultPickle.open("rb") as fh:
            return pickle.load(fh)
    else:
        if showProgress:
            files = list(sourceDir.glob(globPattern))  # Will follow symlinks
            # just like glob(str(sourceDir/globPattern)), but return list[Path]
            for d in tqdm(
                [d for d in sourceDir.iterdir() if d.is_dir()],
                desc=(f'Gathering "{globPattern}" files from {sourceDir.name}'),
            ):
                files.extend(
                    Path(f)
                    for f in glob(
                        str(d / ("**/" + globPattern)), recursive=True
                    )
                )
        else:
            files = [
                Path(f)
                for f in glob(
                    str(sourceDir / ("**/" + globPattern)), recursive=True
                )
            ]
        if resultPickle is not None:
            with resultPickle.open("wb") as fh:
                pickle.dump(files, fh)
    return files
