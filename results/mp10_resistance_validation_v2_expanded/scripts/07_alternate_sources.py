import runpy
from pathlib import Path
ns=runpy.run_path(str(Path(__file__).with_name('01_inputs.py')));globals().update({k:v for k,v in ns.items() if not k.startswith('__')})
jobs=[('Tanaka_supp_actual.pdf','https://insight.jci.org/articles/view/83654/sd/pdf/render/1'),('CLSPN_supp_actual.pdf','https://static-content.springer.com/esm/art%3A10.1007%2Fs00262-023-03388-5/MediaObjects/262_2023_3388_MOESM2_ESM.pdf'),('DTP_TableS6_alt.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/PMC8437243/bin/NIHMS1732990-supplement-14.xlsx?download=1'),('CMP_RNAseq.zip','https://cog.sanger.ac.uk/cmp/download/rnaseq_all_20220624.zip'),('CMP_models.csv','https://cog.sanger.ac.uk/cmp/download/model_list_20240110.csv'),('S1314_code.zip','https://raw.githubusercontent.com/ZB-WCM/GMLF_response_NAC_MIBC/main/GMLF_response_to_NAC_MIBC_npj_Digital_Medicine.zip')]
for item in jobs:get(item)
