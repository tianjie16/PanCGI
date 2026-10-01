import datetime
import json
import subprocess
import sys
import time
from pathlib import Path

root, name, *command = sys.argv[1:]
directory = Path(root)
directory.mkdir(parents=True, exist_ok=True)
record = dict(stage=name, command=command, status='running', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
path = directory/(name+'.json')
if path.exists():
    raise FileExistsError(path)
path.write_text(json.dumps(record, indent=2)+'\n')
start = time.monotonic()
print(f'START {name}', flush=True)
try:
    with (directory/(name+'.log')).open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    record.update(status='pass' if result.returncode == 0 else 'failed', exit_code=result.returncode)
except Exception as exc:
    record.update(status='failed', error=str(exc), exit_code=1)
record.update(elapsed_seconds=time.monotonic()-start, ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
path.write_text(json.dumps(record, indent=2)+'\n')
print(f'{record["status"].upper()} {name}', flush=True)
if record['exit_code']:
    print(f'Failed stage log: {directory/(name+".log")}', file=sys.stderr)
raise SystemExit(record['exit_code'])
