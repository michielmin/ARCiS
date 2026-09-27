#!/usr/bin/env python3
"""Run an ARCiS grid and replace a PyARCiS worker after GGchem fails.

Each worker keeps its initialized opacity tables while it is healthy.  When
GGchem creates ``fatal.case``, the worker reports the failure and exits.  The
supervisor starts a fresh process in the same logical slot and gives the
failed case to that fresh process before assigning it any other work.
"""

from __future__ import annotations

import argparse
import csv
from collections import deque
import math
import multiprocessing as mp
import os
from pathlib import Path
import queue
import shutil
import sys
import time
import traceback
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", default="arcis_noclouds_grid_lumen.dat")
    parser.add_argument("--input", default="output_gas_test/arcis_gas.dat")
    parser.add_argument("--output-root", default="output_gas_test")
    parser.add_argument("--work-root", default="arcis_work")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument(
        "--fatal-retries",
        type=int,
        default=2,
        help="Number of fresh-process retries after the first GGchem failure.",
    )
    parser.add_argument("--planet-mass", type=float, default=1.0)
    parser.add_argument("--threads-per-worker", type=int, default=1)
    parser.add_argument(
        "--home",
        default="/net/lumen/data2/louis",
        help="HOME used by ARCiS to locate its Data directory.",
    )
    parser.add_argument(
        "--run-all",
        action="store_true",
        help="Run rows already marked converged as well.",
    )
    return parser.parse_args()


def read_grid(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        return rows, list(reader.fieldnames)


def write_grid_atomic(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def is_true(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def read_temperature_convergence(path: Path) -> bool:
    if not path.is_file():
        return False
    with path.open() as handle:
        for line in handle:
            key, _, value = line.partition("=")
            if key.strip().lower() == "converged":
                return value.strip().lower() == "true"
    return False


def define_radius(planet_mass_jupiter: float, logg_cgs: float) -> float:
    """Return Rp/Rjup for Mp/Mjup and log10(g [cm s^-2])."""
    gravitational_constant = 6.67430e-8
    jupiter_mass = 1.89816e30
    jupiter_radius = 7.1492e9
    gravity = 10.0**logg_cgs
    return (
        math.sqrt(
            gravitational_constant * planet_mass_jupiter * jupiter_mass / gravity
        )
        / jupiter_radius
    )


def clear_output(output_dir: Path) -> None:
    """Clear generated output but retain the cumulative ARCiS log."""
    for item in output_dir.iterdir():
        if item.name == "log.dat":
            continue
        if item.is_dir() and not item.is_symlink():
            shutil.rmtree(item)
        else:
            item.unlink()


def case_directory(output_root: Path, case: dict[str, str]) -> Path:
    return output_root / (
        f"case_{case['case_id']}"
        f"_Dplanet_{case['Dplanet']}"
        f"_Tint_{case['Tint']}"
        f"_logg_{case['logg']}"
        f"_met_{case['[M/H]']}"
    )


def copy_case_output(
    output_dir: Path,
    output_root: Path,
    case: dict[str, str],
    retry: int,
) -> Path:
    destination = case_directory(output_root, case)
    if destination.exists():
        destination = destination.with_name(
            destination.name + f"_rerun_{retry}_{time.time_ns()}"
        )
    shutil.copytree(output_dir, destination)
    return destination


def flush_and_exit(result_queue: mp.Queue, message: tuple[Any, ...]) -> None:
    """Make sure a fatal/error message reaches the supervisor before exit."""
    result_queue.put(message)
    result_queue.close()
    result_queue.join_thread()


def worker_main(
    token: str,
    task_queue: mp.Queue,
    result_queue: mp.Queue,
    input_file_string: str,
    output_root_string: str,
    work_root_string: str,
    planet_mass: float,
) -> None:
    # Import the extension only inside a freshly spawned process.
    import pyARCiS

    output_root = Path(output_root_string)
    worker_dir = Path(work_root_string) / token
    output_dir = output_root / "worker_output" / token
    worker_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(worker_dir)

    output_dir_string = str(output_dir) + "/"
    try:
        pyARCiS.pyinit(input_file_string, output_dir_string)
        pyARCiS.pyverbose(True)
    except BaseException:
        flush_and_exit(
            result_queue,
            ("init_error", token, traceback.format_exc()),
        )
        return

    result_queue.put(("ready", token))

    while True:
        task = task_queue.get()
        if task is None:
            return

        case, retry = task
        case_id = str(case["case_id"])
        result_queue.put(("started", token, case_id, retry))
        fatal_file = worker_dir / "fatal.case"
        fatal_file.unlink(missing_ok=True)

        try:
            print(
                f"[{token}] case={case_id}, retry={retry}, "
                f"Dplanet={case['Dplanet']}, [M/H]={case['[M/H]']}, "
                f"Tint={case['Tint']}, logg={case['logg']}",
                flush=True,
            )

            pyARCiS.pysetvalue("Mp", float(planet_mass))
            pyARCiS.pysetvalue("Dplanet", float(case["Dplanet"]))
            pyARCiS.pysetvalue("TeffP", float(case["Tint"]))
            pyARCiS.pysetvalue("metallicity", float(case["[M/H]"]))
            radius = define_radius(planet_mass, float(case["logg"]))
            pyARCiS.pysetvalue("Rp", float(radius))

            pyARCiS.pycomputemodel()

            if fatal_file.exists():
                failure_dir = (
                    output_root
                    / "fatal_cases"
                    / f"case_{case_id}_retry_{retry}_{token}"
                )
                failure_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(fatal_file, failure_dir / "fatal.case")
                (failure_dir / "parameters.txt").write_text(
                    "\n".join(f"{key}={value}" for key, value in case.items())
                    + "\n"
                )
                flush_and_exit(
                    result_queue,
                    (
                        "fatal",
                        token,
                        case_id,
                        retry,
                        fatal_file.read_text(errors="replace"),
                    ),
                )
                return

            pyARCiS.pywritefiles()
            destination = copy_case_output(
                output_dir,
                output_root,
                case,
                retry,
            )
            converged = read_temperature_convergence(
                destination / "temperature_convergence.dat"
            )
            clear_output(output_dir)
            result_queue.put(
                ("done", token, case_id, retry, converged, str(destination))
            )
            result_queue.put(("ready", token))

        except BaseException:
            flush_and_exit(
                result_queue,
                ("error", token, case_id, retry, traceback.format_exc()),
            )
            return


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    if args.fatal_retries < 0:
        raise ValueError("--fatal-retries cannot be negative")

    # These values are inherited by spawned workers before PyARCiS is loaded.
    os.environ["HOME"] = args.home
    os.environ["OMP_NUM_THREADS"] = str(args.threads_per_worker)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("PYTHONFAULTHANDLER", "1")

    launch_dir = Path.cwd().resolve()
    grid_file = (launch_dir / args.grid).resolve()
    input_file = (launch_dir / args.input).resolve()
    output_root = (launch_dir / args.output_root).resolve()
    run_id = f"run_{os.getpid()}_{time.time_ns()}"
    work_root = (launch_dir / args.work_root / run_id).resolve()

    if not grid_file.is_file():
        raise FileNotFoundError(grid_file)
    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    rows, fieldnames = read_grid(grid_file)
    required = {"case_id", "Dplanet", "[M/H]", "Tint", "logg"}
    missing = required.difference(fieldnames)
    if missing:
        raise ValueError(f"Grid is missing columns: {sorted(missing)}")
    for optional in ("attempt", "converged"):
        if optional not in fieldnames:
            fieldnames.append(optional)
            for row in rows:
                row[optional] = "0"

    selected = [
        row
        for row in rows
        if args.run_all or not is_true(row.get("converged", "0"))
    ]
    if not selected:
        print("No unconverged cases remain.")
        return 0

    row_by_id = {str(row["case_id"]): row for row in rows}
    if len(row_by_id) != len(rows):
        raise ValueError("case_id values must be unique")

    pending: deque[tuple[dict[str, str], int]] = deque(
        (case, 0) for case in selected
    )
    fresh_retry_by_slot: dict[int, tuple[dict[str, str], int]] = {}
    completed: set[str] = set()
    permanently_failed: set[str] = set()

    context = mp.get_context("spawn")
    result_queue = context.Queue()
    workers: dict[str, dict[str, Any]] = {}
    ready: set[str] = set()
    assigned: dict[str, tuple[dict[str, str], int]] = {}
    inflight: dict[str, tuple[dict[str, str], int]] = {}
    generations = [0 for _ in range(args.workers)]

    def start_worker(slot: int) -> None:
        generation = generations[slot]
        generations[slot] += 1
        token = f"worker_{slot:03d}_generation_{generation:04d}"
        private_queue = context.Queue(maxsize=1)
        process = context.Process(
            target=worker_main,
            args=(
                token,
                private_queue,
                result_queue,
                str(input_file),
                str(output_root),
                str(work_root),
                args.planet_mass,
            ),
            name=token,
        )
        process.start()
        workers[token] = {
            "process": process,
            "queue": private_queue,
            "slot": slot,
        }
        print(f"Started {token} (PID {process.pid})", flush=True)

    def checkpoint(case_id: str, converged: bool) -> None:
        row = row_by_id[case_id]
        row["attempt"] = "1"
        row["converged"] = "1" if converged else "0"
        write_grid_atomic(grid_file, rows, fieldnames)

    def retry_or_fail(
        slot: int,
        task: tuple[dict[str, str], int],
        reason: str,
    ) -> None:
        case, retry = task
        case_id = str(case["case_id"])
        if retry < args.fatal_retries:
            next_task = (case, retry + 1)
            fresh_retry_by_slot[slot] = next_task
            print(
                f"Case {case_id} will be retried by a fresh worker "
                f"({retry + 1}/{args.fatal_retries}): {reason}",
                flush=True,
            )
        else:
            permanently_failed.add(case_id)
            checkpoint(case_id, False)
            print(
                f"Case {case_id} permanently failed after "
                f"{retry + 1} attempts: {reason}",
                flush=True,
            )

    def handle_message(message: tuple[Any, ...]) -> None:
        kind, token, *values = message
        info = workers.get(token)
        if info is None:
            return
        slot = int(info["slot"])

        if kind == "ready":
            ready.add(token)
            return

        if kind == "started":
            task = assigned.pop(token, None)
            if task is not None:
                inflight[token] = task
            return

        if kind == "done":
            case_id, retry, converged, destination = values
            inflight.pop(token, None)
            assigned.pop(token, None)
            completed.add(str(case_id))
            checkpoint(str(case_id), bool(converged))
            print(
                f"Finished case {case_id}; converged={converged}; "
                f"output={destination}",
                flush=True,
            )
            return

        if kind in {"fatal", "error"}:
            case_id, retry, details = values
            task = inflight.pop(token, None) or assigned.pop(token, None)
            ready.discard(token)
            print(
                f"{kind.upper()} in {token}, case {case_id}:\n{details}",
                file=sys.stderr,
                flush=True,
            )
            if task is not None:
                retry_or_fail(slot, task, kind)
            return

        if kind == "init_error":
            (details,) = values
            ready.discard(token)
            print(
                f"Initialization failed in {token}:\n{details}",
                file=sys.stderr,
                flush=True,
            )

    def drain_messages() -> None:
        while True:
            try:
                handle_message(result_queue.get_nowait())
            except queue.Empty:
                return

    def dispatch() -> None:
        for token in list(ready):
            info = workers.get(token)
            if info is None or not info["process"].is_alive():
                ready.discard(token)
                continue
            slot = int(info["slot"])
            if slot in fresh_retry_by_slot:
                task = fresh_retry_by_slot.pop(slot)
            elif pending:
                task = pending.popleft()
            else:
                continue
            info["queue"].put(task)
            assigned[token] = task
            ready.discard(token)

    for slot in range(args.workers):
        start_worker(slot)

    total = len(selected)
    try:
        while len(completed) + len(permanently_failed) < total:
            try:
                handle_message(result_queue.get(timeout=0.5))
            except queue.Empty:
                pass
            drain_messages()

            dead_tokens = [
                token
                for token, info in workers.items()
                if not info["process"].is_alive()
            ]
            for token in dead_tokens:
                workers[token]["process"].join()

            # A child flushes fatal/error messages before exiting.  Drain once
            # more after join so those messages are handled before recovery.
            drain_messages()

            for token in dead_tokens:
                info = workers.pop(token, None)
                if info is None:
                    continue
                slot = int(info["slot"])
                exitcode = info["process"].exitcode
                ready.discard(token)
                task = inflight.pop(token, None) or assigned.pop(token, None)
                info["queue"].close()
                if task is not None:
                    retry_or_fail(
                        slot,
                        task,
                        f"worker exited with code {exitcode}",
                    )

                if len(completed) + len(permanently_failed) < total:
                    start_worker(slot)

            dispatch()

        print(
            f"Grid finished: {len(completed)} completed, "
            f"{len(permanently_failed)} permanently failed.",
            flush=True,
        )
        return 1 if permanently_failed else 0

    finally:
        for token, info in list(workers.items()):
            if info["process"].is_alive():
                try:
                    info["queue"].put(None)
                except (BrokenPipeError, EOFError):
                    pass
        for info in workers.values():
            info["process"].join(timeout=10)
            if info["process"].is_alive():
                info["process"].terminate()
                info["process"].join()


if __name__ == "__main__":
    raise SystemExit(main())
