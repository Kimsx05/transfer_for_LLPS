import runpy
from pathlib import Path
ns=runpy.run_path(str(Path(__file__).with_name('01_inputs.py')));globals().update({k:v for k,v in ns.items() if not k.startswith('__')})
jobs=[('DTP_TableS6_download.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/instance/8437243/bin/NIHMS1732990-supplement-14.xlsx?download=1'),('CrM_publisher.html','https://www.sciencedirect.com/science/article/pii/S0304383524007341'),('CrM_supp1.docx','https://ars.els-cdn.com/content/image/1-s2.0-S0304383524007341-mmc1.docx'),('CrM_supp2.xlsx','https://ars.els-cdn.com/content/image/1-s2.0-S0304383524007341-mmc2.xlsx'),('CMP_models_current.csv','https://cog.sanger.ac.uk/cmp/download/model_list_20240110.csv')]
for item in jobs:get(item)
