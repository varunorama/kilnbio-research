"""Download released ENCODE files and verify their size and published MD5."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import time
import urllib.request
from pathlib import Path


def read_range(url, start, end, size):
    request = urllib.request.Request(url, headers={'Range': f'bytes={start}-{end}'})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                if response.headers.get('Content-Range') != f'bytes {start}-{end}/{size}':
                    raise ValueError('Unexpected range response')
                data = response.read()
            if len(data) != end-start+1:
                raise OSError('Truncated range')
            return data
        except OSError:
            if attempt == 4:
                raise
            time.sleep(2**attempt)


def md5_file(path):
    digest = hashlib.md5()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def fetch(record, destination, workers):
    size = record['file_size']
    md5 = record['md5sum']
    url = record['cloud_metadata']['url']
    final = destination / Path(record['href']).name
    receipt = final.with_suffix(final.suffix+'.md5')
    if final.exists() and final.stat().st_size == size and md5_file(final) == md5:
        receipt.write_text(md5+'\n')
        return
    if receipt.exists():
        receipt.unlink()
    piece = final.with_suffix(final.suffix+'.partial')
    stride = 8 * 1024**2
    def part(start):
        data = read_range(url, start, min(start+stride, size)-1, size)
        with piece.open('r+b') as handle:
            handle.seek(start)
            handle.write(data)
    with piece.open('wb') as handle:
        handle.truncate(size)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(part, range(0,size,stride)))
    if md5_file(piece) != md5:
        raise ValueError(f'ENCODE MD5 mismatch for {final.name}')
    os.replace(piece, final)
    receipt.write_text(md5+'\n')
    print('Verified', final.name, size, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path)
    parser.add_argument('metadata', type=Path)
    parser.add_argument('accessions', nargs='+')
    parser.add_argument('--workers',type=int,default=8)
    args = parser.parse_args()
    args.destination.mkdir(parents=True,exist_ok=True)
    records = {f['accession']:f for p in args.metadata.glob('ENCSR*.json') for f in json.loads(p.read_text())['files']}
    for acc in args.accessions:
        f = records[acc]
        if f['status'] != 'released' or f.get('assembly') != 'GRCh38':
            raise ValueError(acc)
        fetch(f,args.destination,args.workers)
