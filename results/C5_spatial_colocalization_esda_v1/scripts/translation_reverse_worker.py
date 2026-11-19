# Serial worker for predeclared latter half; select available inputs, never effects.
from translation import *
pending=sorted(META.index)[32:][::-1]
while pending:
 available=[s for s in pending if (O/'00_input_audit'/f'{s}.json').exists()]
 if not available:time.sleep(5);continue
 s=available[0];run(s);pending.remove(s)
