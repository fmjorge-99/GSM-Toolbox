"""Background job execution so the UI never freezes.

Long analyses (FBA, FVA, gap-filling, strain design, database builds) are CPU-
bound C code (GLPK/SCIP) that does **not** release Python's GIL, so running them
on a worker *thread* would still block the GUI. Instead each job runs in a
separate **process**; the worker thread simply waits on the result (releasing
the GIL while it waits), keeping the main window fully interactive.

The callable is shipped to the worker process with :mod:`cloudpickle` so the
existing closure/lambda call-sites keep working. If a job cannot be pickled or
the process pool is unavailable, execution falls back to running in-thread so
behaviour is never worse than before.
"""

from __future__ import annotations

import itertools
import json
import multiprocessing as mp
import os
import queue as _queue
import time
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Callable, Dict, List, Optional

from PySide6.QtCore import QObject, QThread, QTimer, Signal

try:
    import cloudpickle
except Exception:  # noqa: BLE001
    cloudpickle = None

# Allow disabling out-of-process execution (e.g. for debugging) via env var.
_PROCESS_MODE = cloudpickle is not None and os.environ.get("GSM_NO_SUBPROCESS") != "1"

# Allow several analyses to run at once (each job is mostly single-core, so cap
# the pool to keep the machine responsive while still using multiple cores).
_MAX_CONCURRENT = max(1, min(4, (os.cpu_count() or 2)))

_executor: ProcessPoolExecutor | None = None


def _get_executor() -> ProcessPoolExecutor | None:
    """Lazily create the reusable worker-process pool (allows concurrent jobs)."""
    global _executor
    if _executor is None:
        try:
            _executor = ProcessPoolExecutor(max_workers=_MAX_CONCURRENT)
        except Exception:  # noqa: BLE001
            return None
    return _executor


def _shutdown_executor() -> None:
    global _executor
    if _executor is not None:
        _executor.shutdown(wait=False, cancel_futures=True)
        _executor = None


def _run_payload(payload: bytes) -> Any:
    """Top-level entry point executed inside the worker process."""
    import cloudpickle as _cp
    fn = _cp.loads(payload)
    return fn()


def _execute(fn: Callable[[], Any]) -> Any:
    """Run ``fn`` in the worker process, falling back to in-thread on any issue."""
    if _PROCESS_MODE:
        try:
            payload = cloudpickle.dumps(fn)
        except Exception:  # noqa: BLE001 - not picklable -> run here
            return fn()
        executor = _get_executor()
        if executor is not None:
            try:
                return executor.submit(_run_payload, payload).result()
            except Exception:  # noqa: BLE001
                # A broken pool can't be reused; drop it and run in-thread.
                _shutdown_executor()
                return fn()
    return fn()


class _Worker(QObject):
    finished = Signal(object)   # result
    failed = Signal(str)        # error message

    def __init__(self, fn: Callable[[], Any]):
        super().__init__()
        self._fn = fn

    def run(self) -> None:
        try:
            result = _execute(self._fn)
        except Exception as exc:  # noqa: BLE001 - surface as a friendly message
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


class JobRunner(QObject):
    """Runs a callable on a background thread (which drives a worker process).

    Keep a reference to the runner alive (e.g. on the main window) until it
    emits ``done``/``error`` — otherwise the thread is garbage-collected.
    """

    done = Signal(object)
    error = Signal(str)

    def __init__(self, parent: QObject | None = None):
        super().__init__(parent)
        self._thread: QThread | None = None
        self._worker: _Worker | None = None

    def start(self, fn: Callable[[], Any]) -> None:
        self._thread = QThread()
        self._worker = _Worker(fn)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._thread.start()

    def _cleanup(self) -> None:
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait()
            self._thread = None
        self._worker = None

    def _on_finished(self, result: Any) -> None:
        self._cleanup()
        self.done.emit(result)

    def _on_failed(self, message: str) -> None:
        self._cleanup()
        self.error.emit(message)


# --------------------------------------------------------------------------- #
# Job manager: tracks several concurrent jobs and estimates their progress.
# --------------------------------------------------------------------------- #
def _durations_path() -> str:
    base = os.path.join(os.path.expanduser("~"), ".gsm_toolbox")
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, "job_durations.json")


# Reasonable starting estimates (seconds) per job kind; refined from real runs.
_DEFAULT_DURATIONS: Dict[str, float] = {
    "fba": 1.0, "pfba": 2.0, "shadow_prices": 1.5, "fva": 8.0,
    "single_reaction_deletion": 20.0, "single_gene_deletion": 25.0,
    "production_envelope": 6.0, "robustness": 6.0, "phase_plane": 12.0,
    "mutant": 4.0, "knockout": 90.0, "overproduction": 10.0, "fseof": 8.0,
    "gimme": 4.0, "eflux": 2.0, "atpm_sensitivity": 5.0, "efm": 15.0,
    "quality_report": 12.0, "blocked_reactions": 10.0,
    "gapfill_growth": 30.0, "gapfill_metabolite": 30.0,
    "download": 30.0, "database": 25.0, "uniprot": 5.0, "generic": 5.0,
}


def _job_entry(payload: bytes, out_queue) -> None:
    """Worker-process entry point: run the cloudpickled callable, return via queue."""
    try:
        import cloudpickle as _cp
        fn = _cp.loads(payload)
        result = fn()
        out_queue.put((True, result))
    except Exception as exc:  # noqa: BLE001
        out_queue.put((False, str(exc)))
    finally:
        try:
            out_queue.close()
            out_queue.join_thread()
        except Exception:  # noqa: BLE001
            pass


class JobManager(QObject):
    """Schedules background jobs across worker processes with a queue.

    Up to ``_MAX_CONCURRENT`` jobs run at once (each in its own process, so they
    use multiple cores and can be individually terminated); additional jobs wait
    in a FIFO **queue** and start automatically as slots free up. Progress is
    *estimated* from how long each job kind has taken before (an EMA persisted
    between sessions), because the solvers expose no step-by-step progress.
    """

    updated = Signal()

    def __init__(self, parent: QObject | None = None):
        super().__init__(parent)
        self._jobs: Dict[int, dict] = {}          # id -> job dict (insertion order = queue order)
        self._ids = itertools.count(1)
        self._durations = self._load_durations()
        self._ctx = mp.get_context("spawn")
        self._timer = QTimer(self)
        self._timer.setInterval(200)
        self._timer.timeout.connect(self._poll)

    # -- persistence ---------------------------------------------------
    def _load_durations(self) -> Dict[str, float]:
        data = dict(_DEFAULT_DURATIONS)
        try:
            with open(_durations_path(), "r", encoding="utf-8") as fh:
                data.update({k: float(v) for k, v in json.load(fh).items()})
        except Exception:  # noqa: BLE001
            pass
        return data

    def _save_durations(self) -> None:
        try:
            with open(_durations_path(), "w", encoding="utf-8") as fh:
                json.dump(self._durations, fh)
        except Exception:  # noqa: BLE001
            pass

    # -- submission ----------------------------------------------------
    def submit(self, fn: Callable[[], Any], *, title: str, kind: str = "generic",
               on_done: Optional[Callable] = None,
               on_error: Optional[Callable] = None) -> int:
        jid = next(self._ids)
        payload = None
        if cloudpickle is not None:
            try:
                payload = cloudpickle.dumps(fn)
            except Exception:  # noqa: BLE001
                payload = None
        self._jobs[jid] = {
            "id": jid, "title": title, "kind": kind, "state": "queued",
            "estimate": float(self._durations.get(kind, self._durations["generic"])),
            "fn": fn, "payload": payload, "on_done": on_done, "on_error": on_error,
            "proc": None, "queue": None, "start": None,
        }
        if not self._timer.isActive():
            self._timer.start()
        self._schedule()
        self.updated.emit()
        return jid

    def _running_count(self) -> int:
        return sum(1 for j in self._jobs.values() if j["state"] == "running")

    def _schedule(self) -> None:
        for job in self._jobs.values():
            if self._running_count() >= _MAX_CONCURRENT:
                break
            if job["state"] == "queued":
                self._start(job)

    def _start(self, job: dict) -> None:
        job["start"] = time.monotonic()
        job["state"] = "running"
        if job["payload"] is None:
            # Cannot ship to a process (cloudpickle missing / unpicklable) — run
            # in-thread as a last resort so the job still completes.
            self._run_inthread(job)
            return
        q = self._ctx.Queue()
        proc = self._ctx.Process(target=_job_entry, args=(job["payload"], q), daemon=False)
        proc.start()
        job["proc"] = proc
        job["queue"] = q

    def _run_inthread(self, job: dict) -> None:
        runner = JobRunner(self)
        job["runner"] = runner
        jid = job["id"]
        runner.done.connect(lambda res, i=jid: self._finish(i, res, success=True))
        runner.error.connect(lambda msg, i=jid: self._finish(i, msg, success=False))
        runner.start(job["fn"])

    # -- polling / completion ------------------------------------------
    def _poll(self) -> None:
        for job in list(self._jobs.values()):
            if job["state"] == "cancelling":
                proc = job.get("proc")
                if proc is None or not proc.is_alive():
                    self._finish(job["id"], "cancelled", success=False, silent=True)
                continue
            if job["state"] != "running" or job.get("queue") is None:
                continue
            try:
                ok, payload = job["queue"].get_nowait()
            except _queue.Empty:
                if job["proc"] is not None and not job["proc"].is_alive():
                    # process ended without a result -> crashed or was terminated
                    if job["state"] == "cancelling":
                        self._finish(job["id"], "cancelled", success=False, silent=True)
                    else:
                        self._finish(job["id"], "The analysis process stopped unexpectedly.",
                                     success=False)
                continue
            self._finish(job["id"], payload, success=ok)
        self._schedule()
        if not self._jobs:
            self._timer.stop()
        self.updated.emit()

    def _finish(self, jid: int, payload, *, success: bool, silent: bool = False) -> None:
        job = self._jobs.pop(jid, None)
        if job is None:
            return
        proc = job.get("proc")
        if proc is not None:
            try:
                if proc.is_alive():
                    proc.terminate()
                proc.join(timeout=2)
            except Exception:  # noqa: BLE001
                pass
        if success and job.get("start") is not None:
            elapsed = time.monotonic() - job["start"]
            prev = self._durations.get(job["kind"], elapsed)
            self._durations[job["kind"]] = max(0.3, 0.7 * prev + 0.3 * elapsed)
            self._save_durations()
        self.updated.emit()
        if silent:
            return
        cb = job["on_done"] if success else job["on_error"]
        if cb is not None:
            cb(payload)

    # -- cancellation ---------------------------------------------------
    def cancel(self, jid: int) -> None:
        job = self._jobs.get(jid)
        if job is None:
            return
        if job["state"] == "queued":
            self._jobs.pop(jid, None)
            self.updated.emit()
            self._schedule()
            return
        # running: terminate its process; the poll loop will reap it silently
        job["state"] = "cancelling"
        proc = job.get("proc")
        if proc is not None:
            try:
                proc.terminate()
            except Exception:  # noqa: BLE001
                pass
        else:
            self._jobs.pop(jid, None)
            self.updated.emit()

    def cancel_all(self) -> None:
        for jid in list(self._jobs.keys()):
            self.cancel(jid)

    # -- queries (for the status bar / dialog) -------------------------
    def count(self) -> int:
        return len(self._jobs)

    def running_count(self) -> int:
        return self._running_count()

    def snapshot(self) -> List[dict]:
        now = time.monotonic()
        out = []
        for job in self._jobs.values():
            running = job["state"] in ("running", "cancelling")
            elapsed = (now - job["start"]) if (running and job["start"]) else 0.0
            est = max(0.5, job["estimate"])
            fraction = min(0.98, elapsed / est) if running else 0.0
            out.append({
                "id": job["id"],
                "title": job["title"],
                "state": job["state"],
                "elapsed": elapsed,
                "estimate": est,
                "fraction": fraction,
                "eta": max(0.0, est - elapsed) if running else est,
            })
        return out
