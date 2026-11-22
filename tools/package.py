"""Package source and evidence without machine-specific compiled binaries."""
from pathlib import Path
import zipfile
import hashlib

root=Path(__file__).resolve().parents[1]
required=[root/'report'/'Group04_report.pdf',root/'slides'/'Group04_slides.pdf',root/'results'/'summary.csv']
for path in required:
    if not path.is_file():
        raise SystemExit(f'Missing deliverable: {path}')
binary_names={'serial','parallel','false_sharing','hello','pingpong'}
archive=root.parent/'Group04_Assignment3.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as output:
    for path in sorted(root.rglob('*')):
        if not path.is_file():
            continue
        if any(part in {'.git','__pycache__'} for part in path.parts):
            continue
        if path.name in binary_names or path.suffix in {'.exe','.o','.pyc'}:
            continue
        output.write(path,Path(root.name)/path.relative_to(root))
digest=hashlib.sha256(archive.read_bytes()).hexdigest()
archive.with_suffix('.zip.sha256').write_text(f'{digest}  {archive.name}\n')
print(f'{archive}: {archive.stat().st_size:,} bytes; SHA256 {digest}')
