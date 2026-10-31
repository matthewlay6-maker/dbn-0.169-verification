"""Cloud Run Job entrypoint: run this task's work items with /app/dbn, streaming every output line
to stdout (Cloud Logging) and to /out/task-NN.jsonl (GCS volume, copied at least once a minute)."""
import json, os, subprocess, shutil, time

idx = int(os.environ.get("CLOUD_RUN_TASK_INDEX", "0"))
cfg = json.load(open("/app/tasks.json"))
prefix = os.environ.get("OUT_PREFIX", "/out")
os.makedirs(prefix, exist_ok=True)
out_local = f"/tmp/task-{idx:03d}.jsonl"
out_remote = f"{prefix}/task-{idx:03d}.jsonl"
last_copy = 0.0


def emit(line, force=False):
    global last_copy
    line = line.rstrip("\n")
    print(f"[task {idx}] {line}", flush=True)
    with open(out_local, "a") as fh:
        fh.write(line + "\n")
    if force or time.time() - last_copy > 60:
        shutil.copyfile(out_local, out_remote)
        last_copy = time.time()


emit(json.dumps({"task": idx, "start": time.time(), "items": cfg["work"][idx]}), force=True)
for item in cfg["work"][idx]:
    if item[0] in ("block", "blockpart"):
        cmd = ["/app/dbn", "block", str(item[1]), str(item[2]), cfg["t0"], cfg["y0"]] + ([str(item[3]), str(item[4])] if item[0] == "blockpart" else [])
    else:
        cmd = ["/app/dbn", "barrier", cfg["X0"], cfg["y0"], cfg["t0"], repr(item[1]), repr(item[2]),
               str(cfg["level"]), cfg["margin"]]
    t0 = time.time()
    env = dict(os.environ)
    if len(item) > 5 and isinstance(item[5], dict):
        env.update({k: str(v) for k, v in item[5].items()})
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    for line in p.stdout:
        emit(line)
    rc = p.wait()
    emit(json.dumps({"item": item, "returncode": rc, "seconds": time.time() - t0,
                     "stderr_tail": p.stderr.read()[-500:]}), force=True)
emit(json.dumps({"task": idx, "done": True, "end": time.time()}), force=True)
