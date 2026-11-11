from pathlib import Path
import csv,re
p=Path('/home/jinshengxi/bladder_metabolism/results/llps_C5_cellchat_v2')
for f in sorted((p/'00_audit').glob('worker_*.log')):
 s=f.read_text(errors='replace');st=s[s.rfind('START '):];q=re.findall(r'\|\s*(\d+)%',st)
 print(f.name,st.splitlines()[0] if st else '',q[-1]+'%' if q else 'preparing')
print('Unfiltered inference caches:',len(list((p/'00_audit/unfiltered').glob('*.rds'))),'/59')
for d in sorted(p.glob('0[123]_threshold*')):print(d.name,'complete',len(list(d.glob('*/SUCCESS.csv'))),'/59')
