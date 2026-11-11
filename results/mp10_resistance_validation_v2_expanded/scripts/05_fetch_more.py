import runpy
from pathlib import Path
ns=runpy.run_path(str(Path(__file__).with_name('01_inputs.py')));globals().update({k:v for k,v in ns.items() if not k.startswith('__')})
jobs=[('DTP.xml','https://pmc-oa-opendata.s3.amazonaws.com/PMC8437243.1/PMC8437243.1.xml'),('Tanaka_supp.pdf','https://insight.jci.org/articles/view/83654/sd/1'),('GDSC_expression.zip','https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/Cell_line_RMA_proc_basalExp.txt.zip'),('cissig_tree.json','https://api.github.com/repos/jessicascarborough/cissig/git/trees/c8317dc0a6ab6b8ee66dc98f2e817014f1a3c2a5?recursive=1'),('CrM_article.html','https://www.sciencedirect.com/science/article/pii/S0304383524007328')]
for item in jobs:get(item)
