import os,resource,sys,runpy,time,json
from pathlib import Path
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
resource.setrlimit(resource.RLIMIT_AS,(100_000_000_000,100_000_000_000))
p=sys.argv.pop(1);sys.path.insert(0,os.path.dirname(os.path.abspath(p)));sys.argv[0]=p;t=time.time()
try:runpy.run_path(p,run_name='__main__')
finally:
 r=resource.getrusage(resource.RUSAGE_SELF)
 out=Path(__file__).resolve().parents[1]/'logs'/('resources_'+Path(p).stem+'.json')
 out.write_text(json.dumps(dict(script=p,elapsed_seconds=time.time()-t,max_rss_kb=r.ru_maxrss,user_cpu_seconds=r.ru_utime,system_cpu_seconds=r.ru_stime,cpu_affinity=sorted(os.sched_getaffinity(0)),threads=4,memory_limit_bytes=100_000_000_000),indent=2))
