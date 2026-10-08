"""Download verified official human NicheNet v2 resources without inference."""
import concurrent.futures
import hashlib
import json
import pathlib
import tarfile
import urllib.request
import urllib.parse

ROOT = pathlib.Path('/home/jinshengxi/bladder_metabolism/cellchat_nichechat_fibro7')
OUT = ROOT / '05_nichenet_preflight'
OUT.mkdir(exist_ok=True)
HEADERS = {'User-Agent': 'fibro7-nichenet-preflight/1.0'}

def read_url(url):
    request = urllib.request.Request(url, headers=HEADERS)
    if urllib.parse.urlparse(url).hostname == 'zenodo.org':
        # Existing local proxy fails the CERN TLS connection; use direct TLS
        # only for the verified official resource domain. Certificates verified.
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        return opener.open(request, timeout=120)
    return urllib.request.urlopen(request, timeout=120)

def download(url, path, expected_md5=None):
    if path.exists():
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        if expected_md5 and digest != expected_md5:
            raise RuntimeError(f'Existing file checksum mismatch: {path}')
        print(f'REUSED {path.name}', flush=True)
        return
    part = path.with_name(path.name + '.part')
    with read_url(url) as response, part.open('wb') as handle:
        digest = hashlib.md5()
        total = 0
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            handle.write(block)
            digest.update(block)
            total += len(block)
    if expected_md5 and digest.hexdigest() != expected_md5:
        raise RuntimeError(f'Download checksum mismatch: {path}')
    part.rename(path)
    print(f'DOWNLOADED {path.name}: {total} bytes; MD5={digest.hexdigest()}', flush=True)

commit_path = OUT / 'official_nichenetr_commit.json'
if commit_path.exists():
    commit = json.loads(commit_path.read_text())
else:
    with read_url('https://api.github.com/repos/saeyslab/nichenetr/commits/master') as handle:
        commit = json.load(handle)
    commit_path.write_text(json.dumps(commit, indent=2))
sha = commit['sha']
source = OUT / f'nichenetr_{sha}.tar.gz'
download(f'https://codeload.github.com/saeyslab/nichenetr/tar.gz/{sha}', source)
if not (OUT / f'nichenetr-{sha}').exists():
    with tarfile.open(source) as archive:
        archive.extractall(OUT, filter='data')
print(f'OFFICIAL SOURCE READY commit={sha}', flush=True)

with read_url('https://zenodo.org/api/records/7074291') as handle:
    record = json.load(handle)
(OUT / 'official_zenodo_7074291_metadata.json').write_text(json.dumps(record, indent=2))
required = {'ligand_target_matrix_nsga2r_final.rds', 'lr_network_human_21122021.rds', 'weighted_networks_nsga2r_final.rds'}
files = [item for item in record['files'] if item['key'] in required]
if {item['key'] for item in files} != required:
    raise RuntimeError('Official record does not contain all three documented human resources')
manifest = []
for item in files:
    checksum = item['checksum']
    if not checksum.startswith('md5:'):
        raise RuntimeError('Unexpected official checksum algorithm')
    entry = {'resource': item['key'], 'url': item['links']['self'], 'official_checksum': checksum, 'size_bytes': item['size'], 'local_path': str(OUT / item['key']), 'record': '7074291', 'source_commit': sha}
    manifest.append(entry)
    print('RESOURCE', entry, flush=True)
(OUT / 'resource_download_manifest.json').write_text(json.dumps(manifest, indent=2))
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    futures = [pool.submit(download, entry['url'], pathlib.Path(entry['local_path']), entry['official_checksum'].split(':', 1)[1]) for entry in manifest]
    for future in futures:
        future.result()
print('ALL OFFICIAL NICHENET RESOURCES AVAILABLE', flush=True)
