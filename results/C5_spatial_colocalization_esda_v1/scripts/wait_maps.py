import time,subprocess
from pathlib import Path
p=Path('results/C5_spatial_colocalization_esda_v1')
while len(list((p/'03_local_bivariate_moran').glob('*.parquet')))<63:time.sleep(5)
subprocess.run(['bash',str(p/'scripts/run_R.sh'),str(p/'scripts/plots.R'),'maps'],check=True)
