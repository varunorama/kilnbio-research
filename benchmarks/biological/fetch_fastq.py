"""Fetch an ENA-derived TSV (url, file, bytes, md5); resume and verify every file.

Uses curl for resumable HTTPS transfers. Outputs and manifests belong outside Git.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import http.client
import json
from pathlib import Path
import re
import subprocess
import time
import urllib.request
from urllib.parse import urlsplit


def digest(path):
    h = hashlib.md5()
    with path.open('rb') as src:
        for block in iter(lambda: src.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def ranged_transfer(url, part, size, workers):
    """Persist completed 8 MiB ranges; only a final MD5 permits publication."""
    stride=8*1024**2
    chunks=part.with_name(part.name+'.chunks')
    chunks.mkdir(exist_ok=True)
    # Reuse the complete prefix from an interrupted sequential transfer.
    if part.exists() and not (chunks/'assembled').exists():
        with part.open('rb') as src:
            for start in range(0,min(part.stat().st_size,size),stride):
                data=src.read(min(stride,size-start))
                target=chunks/str(start)
                if len(data)==min(stride,size-start) and not target.exists():
                    target.write_bytes(data)
    def get(start):
        end=min(start+stride,size)-1;target=chunks/str(start)
        if target.exists() and target.stat().st_size==end-start+1:return
        for attempt in range(6):
            try:
                req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'})
                with urllib.request.urlopen(req,timeout=90) as resp:
                    if resp.headers.get('Content-Range')!=f'bytes {start}-{end}/{size}':
                        raise ValueError('Server returned an incorrect byte range')
                    data=resp.read()
                if len(data)!=end-start+1: raise OSError('Truncated range')
                tmp=target.with_suffix('.tmp');tmp.write_bytes(data);tmp.replace(target)
                return
            except (OSError,ValueError,http.client.HTTPException):
                if attempt==5:raise
                time.sleep(2**attempt)
    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        list(pool.map(get,range(0,size,stride)))
    (chunks/'assembled').touch()
    with part.open('wb') as dst:
        for start in range(0,size,stride):dst.write((chunks/str(start)).read_bytes())
    # Retain chunks until the full-file checksum has passed.
    return chunks


def fetch(row, out, range_workers=1):
    name, size, checksum = row['file'], int(row['bytes']), row['md5']
    if Path(name).name != name or not re.fullmatch(r'[A-Za-z0-9_.-]+', name):
        raise ValueError('Unsafe destination name')
    parsed = urlsplit(row['url'])
    if parsed.scheme != 'https' or parsed.hostname != 'ftp.sra.ebi.ac.uk' or parsed.username:
        raise ValueError('Expected public ENA HTTPS URL')
    if size <= 0 or not re.fullmatch('[0-9a-f]{32}', checksum):
        raise ValueError('Missing ENA size/MD5')
    final, part = out / name, out / (name + '.part')
    if final.exists():
        if final.stat().st_size != size or digest(final) != checksum:
            raise ValueError(f'Existing file failed integrity: {final}')
    else:
        print('Downloading', name, size, flush=True)
        command = ['curl', '--fail', '--location', '--silent', '--show-error',
                   '--retry', '8', '--retry-delay', '5',
                   '--connect-timeout', '60', '--speed-limit', '1024', '--speed-time', '180',
                   '--continue-at', '-', '--output', str(part), row['url']]
        chunks=ranged_transfer(row['url'],part,size,range_workers) if range_workers>1 else None
        if chunks is None:subprocess.run(command, check=True)
        if part.stat().st_size != size or digest(part) != checksum:
            raise ValueError(f'Download failed integrity: {part}')
        part.replace(final)
        if chunks:
            for start in range(0,size,8*1024**2):
                (chunks/str(start)).unlink()
            (chunks/'assembled').unlink()
            if not any(chunks.iterdir()):chunks.rmdir()
    receipt = {'url': row['url'], 'bytes': size, 'md5': checksum, 'verified': True}
    final.with_name(name + '.verified.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print('Verified', name, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest', type=Path)
    p.add_argument('destination', type=Path)
    p.add_argument('--workers', type=int, choices=range(1, 9), default=3)
    p.add_argument('--range-workers', type=int, choices=range(1, 17), default=1)
    args = p.parse_args()
    with args.manifest.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows or len({r['file'] for r in rows}) != len(rows):
        raise ValueError('Manifest empty or duplicated')
    args.destination.mkdir(parents=True, exist_ok=True)
    failures=[]
    with concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        pending={pool.submit(fetch,r,args.destination,args.range_workers):r['file'] for r in rows}
        for task in concurrent.futures.as_completed(pending):
            try:task.result()
            except Exception as e:
                failures.append(f'{pending[task]}: {e}')
    if failures:p.exit(1,'Failed downloads; verified files and completed ranges are retained:\n'+'\n'.join(failures)+'\n')


if __name__ == '__main__':
    main()
