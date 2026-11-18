from pathlib import Path
import time,subprocess
O=Path('results/C5_fivegene_LLPS77_spatial_v1');ready=O/'tables/matched_spot_pairs.parquet'
while not ready.exists():time.sleep(5)
time.sleep(5)
for name in ['threshold_rings.py','summarize.py','supplement.py']:
 print('START',name,flush=True);subprocess.run(['bash',str(O/'scripts/run_python.sh'),str(O/'scripts'/name)],check=True)
(O/'00_input_audit/POST_STATISTICS_COMPLETE.txt').write_text('Primary and threshold sensitivity, patient/dataset summaries and additional controls complete.\n')
print('POST STATISTICS COMPLETE',flush=True)
