"""Interruptible local supervisor. Heavy child process group runs only 23:00-07:00 CST."""
import os,sys,time,signal,subprocess,json,fcntl
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))  # shared hard cap inherited by all descendants
r=Path(__file__).resolve().parents[1];lock=open(r/'logs/night_worker.lock','w')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:sys.exit('A night supervisor is already active')
child=None;paused=False;last=None
try:
 while True:
  now=datetime.now(ZoneInfo('Asia/Shanghai'));allowed=os.environ.get('MP10_RUN_NOW')=='1' or now.hour>=23 or now.hour<7
  if child is None and allowed:child=subprocess.Popen(sys.argv[1:],start_new_session=True);print('HEAVY_STARTED',now.isoformat(),child.pid,flush=True)
  if child is not None:
   if child.poll() is not None:
    (r/'logs/night_status.json').write_text(json.dumps(dict(time=now.isoformat(),state='worker_complete' if child.returncode==0 else 'worker_finished_with_errors',returncode=child.returncode,child_pid=child.pid,supervisor_pid=os.getpid()),indent=2))
    print('HEAVY_EXIT',now.isoformat(),child.returncode,flush=True);sys.exit(child.returncode)
   if not allowed and not paused:os.killpg(child.pid,signal.SIGSTOP);paused=True;print('PAUSED_AT_WINDOW_END',now.isoformat(),flush=True)
   elif allowed and paused:os.killpg(child.pid,signal.SIGCONT);paused=False;print('RESUMED_IN_WINDOW',now.isoformat(),flush=True)
  state='waiting_for_23_CST' if child is None else ('paused_until_23_CST' if paused else 'running')
  (r/'logs/night_status.json').write_text(json.dumps(dict(time=now.isoformat(),state=state,child_pid=child.pid if child else None,supervisor_pid=os.getpid()),indent=2))
  if state!=last:print(state,now.isoformat(),flush=True);last=state
  time.sleep(20)
except KeyboardInterrupt:
 if child is not None:
  os.killpg(child.pid,signal.SIGCONT);os.killpg(child.pid,signal.SIGTERM)
 raise
