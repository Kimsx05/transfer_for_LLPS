from pathlib import Path
import json
import numpy as np,pandas as pd,pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path('results/C5_spatial_exploratory_all63_v1');cols=['major_mean__Fibroblast','major_mean__Epithelial','MMP11_log1p_CP10k','nCount'];v={k:[] for k in cols}
files=sorted((O/'01_spatial_scores').glob('*.parquet'));assert len(files)==63
for f in files:
 d=pd.read_parquet(f,columns=cols)
 for k in cols:v[k].append(d[k].to_numpy())
(O/'00_inputs_and_QC/shared_color_limits.json').write_text(json.dumps({k:float(np.quantile(np.concatenate(x),.99)) for k,x in v.items()},indent=2))
