# For RNA-Seq alignment, with each file or pair of files per sample. Out put aligned file per sample.


from subprocess import run
from tempfile import NamedTemporaryFile
from Bio import SeqIO
import os
import time
import argparse
import logging



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw', required=True, nargs='+', help='path(s) to raw data')
    parser.add_argument('--out', required=True, help='path to output alignment')
    parser.add_argument('--genome', nargs="+", required=True, help='path to genome file(s), also supports indexed genome (*.bt2)')

    parser.add_argument('--isPe', action='store_true', help='set if you have pairend')

    parser.add_argument('--rawExt', help='extension ".fastq.gz" or ".fq.gz"')
    parser.add_argument('--pesuffix', nargs=2, help='suffix to pairend file name. eg. a_1_.fq.gz and a_2_.fq.gz, then set this to --pesuffix _1_ _2_')
    parser.add_argument('--ncpu', type=str, help='number of cpu to use')
    parser.add_argument('--sampleNames', nargs="+", help='sample names if file names contain random generated string which will unpair the samples after removing pairend suffix')

    args = parser.parse_args()

    # required arguments
    rawFolders = args.raw
    outFolder = args.out
    genomes = args.genome

    # boolean
    isPe = args.isPe
    peSfx = args.pesuffix

    # optional
    rawExt = args.rawExt
    samples = args.sampleNames
    ncpu = args.ncpu

    if isinstance(ncpu, type(None)):
        ncpu = "1" 

    #print(args)
    #exit()

    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
    logging.basicConfig(filename=os.path.join(outFolder, 'align.log'), level=logging.DEBUG)


    toBamNcpu = str(max(int(int(ncpu)/8), 1))


    logging.debug(args)

    logging.info('='*20 + getTime() + '='*20)

    genomeBowtie2Idx = buildBowtie2idx(genomes, out=os.path.join(outFolder, 'genomeIdx'))

    # Gether files for each sample
    files = []
    for folder in rawFolders:
        for f in os.listdir(folder):
            # rawExt defined here if not defined
            if isinstance(rawExt, type(None)):
                fn, ext = os.path.splitext(f)
                if ext == '.gz':
                    ext = os.path.splitext(fn)[1] + ext
                elif ext in ['.fastq', '.fq']:
                    pass
                else:
                    continue
                rawExt = ext
                logging.info(f'Found extension for raw files {rawExt}')
            if f.endswith(rawExt):
                files.append(os.path.join(folder,f))
    assert len(files) > 0, f'Files not found in {rawFolders}'
    files.sort()
            
    # peSfx defined here if not defined
    if isPe:
        peSfx = imputePesuffix(files, rawExt, peSfx)

    # samples defined here if not defined
    if isinstance(samples, type(None)):
        fns = [os.path.split(f)[1][:-len(rawExt)] for f in files]
        if isPe:
            samples = list(set(f[:-len(peSfx[0])] for f in fns))
        else:
            samples = fns
        samples.sort()


    sampleFileDict = {}
    for s in samples:
        fs = [f for f in files if s in f]
        fs.sort()
        if isPe:
            sampleFileDict[s] = {peSfx[0]: [f for f in fs if f[:-len(rawExt)].endswith(peSfx[0])],
                                 peSfx[1]: [f for f in fs if f[:-len(rawExt)].endswith(peSfx[1])]}
        else:
            sampleFileDict[s] = fs 

    totalTs = time.time()
    assert len(samples) > 0
    logging.info(f'Samples to process: {samples}')
    for i, s in enumerate(samples):
        ts = time.time()
        logging.info(f'Processing {i+1}/{len(samples)}: {s}')

        # prepare align arguments
        args = ['bowtie2', '-x', genomeBowtie2Idx, '-p', ncpu]
        if isPe:
            args.extend([
                '-1', ','.join(sampleFileDict[s][peSfx[0]]),
                '-2', ','.join(sampleFileDict[s][peSfx[1]])
            ])
        else:
            args.extend([
                '-U', ','.join(sampleFileDict[s])
            ])

        # prepare convert to bam arguments
        target = os.path.join(outFolder, f'{s.strip("_")}.bam')

        if os.path.isfile(target):
            pass
        args.extend([
            '|', 'samtools', 'view', '-bS', '-@', toBamNcpu,
            "|", "samtools", "sort", '-@', toBamNcpu, "--write-index", '-o', target
        ])
        logging.info(' '.join(args))

        # Start running both
        result = run(args=' '.join(args), shell=True, capture_output=True)
        if result.returncode != 0:
            logging.info('stderr: ' + result.stderr.decode())
            logging.info('stdout: ' + result.stdout.decode())
        logging.info(result.stderr.decode()) # stderr has logging.info info from bowtie2
        logging.info(f'Finished in {diffTime(ts)}\n')

    logging.info(f'All done, time elapsed {diffTime(totalTs)}')
    logging.info('='*20 + getTime() + '='*20 + '\n'*2)


def buildBowtie2idx(fs, out, name=None):
    if type(fs) == str:
        fs = [fs]
    fs = [os.path.realpath(f) for f in fs]

    # check if idx file exists
    genomePath, genomeName = os.path.split(fs[0])
    genomeName = os.path.splitext(genomeName)[0]
    logging.info(f'genome name: {genomeName}')

    out = os.path.realpath(out)
    if not os.path.isdir(out):
        os.makedirs(out)
    bt2_base = ('_'.join(os.path.splitext(os.path.split(f)[1])[0] for f in fs) if name == None else name)
    outIdxForUse = os.path.join(out, bt2_base)

    # see if index file already exists
    outIdxFile0 = outIdxForUse + '.1.bt2'
    outIdxLargeFile0 = outIdxForUse + '.1.bt21'
    if os.path.isfile(outIdxFile0):
        logging.info(f'Found {outIdxFile0}, using it as index.')
        return outIdxForUse
    elif os.path.isfile(outIdxLargeFile0):
        logging.info(f'Found {outIdxLargeFile0}, using it as index.')
        return outIdxForUse
    
    # convert gbk to fa
    converted = []
    if os.path.splitext(fs[0])[1] not in ['.fa', '.fasta', '.fna', '.fsa']:
        # try gbk
        newFs = []
        try:
            for f in fs:
                for s in SeqIO.parse(f, 'genbank'):
                    newF = NamedTemporaryFile()
                    SeqIO.write(s, newF.name, 'fasta')
                    newFs.append(newF.name)
                    converted.append(newF) # for closing these files later
        except BaseException as err:
            logging.error(f'Unexpected error {err=}, {type(err)=}') # python >= 3.8
        fs = newFs

    logging.info('-'*20 + 'Indexing genome ' + getTime() + '-'*20)
    args = [
        'bowtie2-build',
        ','.join(fs),
        os.path.join(out,bt2_base),
    ]
    logging.info(' '.join(args))
    result = run(args=args, capture_output=True)
    for f in converted:
        f.close()
    if result.returncode != 0:
        logging.info('stderr: ' + result.stderr.decode())
        logging.info('stdout: ' + result.stdout.decode())
        logging.info('-'*20 + 'Error Indexing genome' + getTime() + '-'*20)
        raise Exception
    logging.info('-'*20 + 'DONE Indexing genome' + getTime() + '-'*20)
    logging.info('\n'*2)

    return os.path.join(out, bt2_base)


def imputePesuffix(rawFiles, rawExt, peSfx=None):
    if not isinstance(peSfx, type(None)):
        return peSfx
    assert len(rawFiles)%2 == 0
    assert len(rawFiles) > 0
    fns = [f[:-len(rawExt)] for f in rawFiles]
    fns.sort()
    for i in range(1, 12):
        s = set(fn[-i:] for fn in fns)
        if len(s) == 2:
            peSfx = list(s)
            peSfx.sort()
            logging.info(f'Found pairend suffix {peSfx}')
            return peSfx
    raise ValueError(f'pair end suffix not found in {rawFiles}')





def getTime():
    return time.strftime('%z, %a, %d %b %Y, %H:%M:%S', time.localtime())


def diffTime(a):
    d = abs(time.time()-a)
    h = int(d//3600)
    return str(h).zfill(2) + time.strftime(':%M:%S', time.gmtime(d))



if __name__ == "__main__":
    main()
