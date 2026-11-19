"""Install a user-local KilnBio environment and verify its first analysis."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MAMBA_VERSION = "2.10.0-0"
MAMBA_SHA256 = "701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd"
PROFILES = {
    "core": ["python=3.12", "pip"],
    "bladder": ["r-base=4.6.2", "c-compiler", "cxx-compiler", "fortran-compiler",
                "make", "pkg-config", "libcurl", "libxml2", "libxml2-devel", "openssl", "zlib", "bzip2", "xz"],
    "analysis": [],
    "singlecell": [],
    "sequencing": [],
    "assistant": [],
}
REQUIREMENTS = {
    "core": "environments/report-requirements.txt",
    "analysis": "environments/multiomic-requirements.txt",
    "singlecell": "environments/analysis-requirements.txt",
    "assistant": "assistants/requirements.txt",
}


def profiles_requested(names):
    names = set(names or ["core"])
    if "all" in names:
        names = set(PROFILES)
    names.add("core")
    if "singlecell" in names:
        names.add("analysis")
    if "analysis" in names:
        names.add("bladder")
    return sorted(names)


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def source_digest():
    files = sorted(p for p in ROOT.rglob("*") if p.is_file()
                   and not any(part in (".git", "__pycache__", ".venv") for part in p.relative_to(ROOT).parts)
                   and p.suffix != ".pyc")
    return hashlib.sha256("\n".join(str(p.relative_to(ROOT)) + ":" + digest(p)
                                    for p in files).encode()).hexdigest()


def execute(args, **kwargs):
    print("Running " + shlex.join([str(x) for x in args]), flush=True)
    return subprocess.run([str(x) for x in args], check=True, **kwargs)


def container_status(engine):
    binary = shutil.which(engine)
    if not binary:
        return "missing"
    args = [binary, "info"] if engine == "docker" else [binary, "--version"]
    try:
        result = subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                timeout=30)
        return "ready" if result.returncode == 0 else "unavailable"
    except (OSError, subprocess.TimeoutExpired):
        return "unavailable"


def smoke(launcher):
    execute([launcher, "--help"], stdout=subprocess.DEVNULL)
    with tempfile.TemporaryDirectory(prefix="kilnbio-install-check-") as directory:
        folder = Path(directory)
        roster = folder / "roster.tsv"
        shutil.copy2(ROOT / "examples/assay_roster.tsv", roster)
        result = folder / "reach.json"
        execute([launcher, "assay-audit", "reach", roster, "--out", result])
        values = json.loads(result.read_text())
        if values.get("donors") != 4 or values.get("base_reach") != 3 or values.get("combined_reach") != 2:
            raise RuntimeError("The installed donor audit returned unexpected results")
    return {"donors": 4, "base_reach": 3, "combined_reach": 2}


def check(prefix):
    receipt = json.loads((prefix / "installation.json").read_text())
    smoke(prefix / "bin/kilnbio")
    if "sequencing" in receipt["profiles"]:
        status = container_status(receipt["container"])
        print("Container engine " + receipt["container"] + ": " + status)
        if status != "ready":
            raise RuntimeError("Enable the selected container engine before running FASTQs. See docs/INSTALLATION.md")
    print("Installation checks passed")


def install(args, profiles):
    prefix = args.prefix.expanduser().resolve()
    if prefix == ROOT or ROOT in prefix.parents:
        raise ValueError("Choose an installation prefix outside the source checkout")
    prefix.mkdir(parents=True, exist_ok=True)
    lock = prefix / ".setup-lock"
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError("Another setup may be running. If it has stopped, remove " + str(lock))
    try:
        install_locked(args, profiles, prefix)
    finally:
        lock.rmdir()


def install_locked(args, profiles, prefix):
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
        raise ValueError("Automatic environments currently support Linux x86-64 and WSL2. See docs/INSTALLATION.md for other systems")
    if not shutil.which("curl") or not shutil.which("bash"):
        raise ValueError("Install curl and bash before setup")
    packages = sorted({package for name in profiles for package in PROFILES[name]})
    requirement_files = [ROOT / REQUIREMENTS[name] for name in profiles if name in REQUIREMENTS]
    for file in requirement_files:
        if not file.is_file():
            raise ValueError("This checkout lacks " + str(file.relative_to(ROOT)))
    source_hash = source_digest()
    spec = {"profiles": profiles, "packages": packages, "source_sha256": source_hash}
    identifier = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]
    bundle = prefix / "bundles" / identifier
    bundle.mkdir(parents=True, exist_ok=True)
    tools = prefix / "tools" / MAMBA_VERSION
    tools.mkdir(parents=True, exist_ok=True)
    mamba = tools / "micromamba"
    if not mamba.is_file() or digest(mamba) != MAMBA_SHA256:
        partial = tools / (mamba.name + ".download")
        url = "https://github.com/mamba-org/micromamba-releases/releases/download/" + MAMBA_VERSION + "/micromamba-linux-64"
        execute(["curl", "--fail", "--location", "--retry", "3", "--connect-timeout", "20",
                 "--max-time", "300", url, "-o", partial])
        if digest(partial) != MAMBA_SHA256:
            raise ValueError("Micromamba checksum mismatch")
        partial.chmod(0o755)
        partial.replace(mamba)
    environment_spec = {"profiles": profiles, "packages": packages,
                        "r_setup": digest(ROOT / "environments/setup_r.R"),
                        "bladder_setup": {name: digest(ROOT / "recipes/bladder" / name) for name in ("setup.py", "setup.R")},
                        "requirements": {str(p.relative_to(ROOT)): digest(p) for p in requirement_files}}
    environment_id = hashlib.sha256(json.dumps(environment_spec, sort_keys=True).encode()).hexdigest()[:16]
    environment = prefix / "environments" / environment_id
    mamba_env = os.environ.copy()
    mamba_env["MAMBA_ROOT_PREFIX"] = str(prefix / "mamba")
    r_library = environment / "kilnbio-R"
    mamba_env["R_LIBS_USER"] = str(r_library)
    base = [mamba, "--no-rc"]
    run = base + ["run", "-p", environment]
    ready = prefix / "environments" / (environment_id + ".json")
    if not ready.exists():
        action = "install" if (environment / "conda-meta").is_dir() else "create"
        execute(base + [action, "-y", "-p", environment, "--override-channels", "-c", "conda-forge",
                        "-c", "bioconda", "--strict-channel-priority", *packages], env=mamba_env)
        pip_args = [value for file in requirement_files for value in ("-r", str(file))]
        if pip_args:
            execute(run + ["python", "-m", "pip", "install", *pip_args], env=mamba_env)
        execute(run + ["python", "-m", "pip", "check"], env=mamba_env)
        if "bladder" in profiles:
            if "analysis" in profiles:
                execute(run + ["Rscript", "--vanilla", ROOT / "environments/setup_r.R", r_library], env=mamba_env)
            pin_script = ROOT / "recipes/bladder/setup.py"
            if not pin_script.is_file():
                raise ValueError("This checkout needs recipes/bladder/setup.py for the pinned R packages")
            pinned = environment / "bladder-R"
            if not (pinned / "kilnbio-pinned.json").is_file():
                # A fresh staging library allows a failed source build to be retried safely.
                with tempfile.TemporaryDirectory(prefix="bladder-r-", dir=bundle) as staging:
                    library = Path(staging) / "library"
                    execute(run + ["python", pin_script, "--library", library], env=mamba_env)
                    if pinned.exists():
                        raise ValueError("Unexpected incomplete pinned library at " + str(pinned))
                    library.replace(pinned)
                    (pinned / "kilnbio-pinned.json").write_text(json.dumps({"limma": "3.68.5", "statmod": "1.5.2"}))
            mamba_env["R_LIBS_USER"] = str(pinned) + os.pathsep + str(r_library)
            execute(run + ["Rscript", "--vanilla", "-e",
                    'stopifnot(getRversion()=="4.6.2", packageVersion("limma")=="3.68.5", packageVersion("statmod")=="1.5.2")'], env=mamba_env)
        ready.write_text(json.dumps(environment_spec, indent=2) + "\n")
    with (bundle / "conda-explicit.txt").open("w") as stream:
        execute(base + ["list", "-p", environment, "--explicit"], env=mamba_env, stdout=stream)
    with (bundle / "pip-freeze.txt").open("w") as stream:
        execute(run + ["python", "-m", "pip", "freeze"], env=mamba_env, stdout=stream)
    code_prefix = bundle / "code"
    execute(run + ["bash", ROOT / "install.sh", code_prefix, "no-runtime"], env=mamba_env,
            stdout=subprocess.DEVNULL)
    if "sequencing" in profiles:
        runtime = bundle / "runtime"
        if not (runtime / "kilnbio-ready.json").exists():
            execute(run + ["bash", ROOT / "setup_linux_runtime.sh", runtime], env=mamba_env)
            (runtime / "kilnbio-ready.json").write_text(json.dumps({"java": 21, "nextflow": "25.10.4"}))
    wrapper = bundle / "kilnbio"
    lines = ["#!/usr/bin/env bash", "set -Eeuo pipefail",
             "export MAMBA_ROOT_PREFIX=" + shlex.quote(str(prefix / "mamba")),
             "export KILNBIO_MANAGED_LAUNCHER=" + shlex.quote(str(wrapper)),
             "export KILNBIO_INSTALL_PREFIX=" + shlex.quote(str(prefix))]
    if "sequencing" in profiles:
        lines.append("source " + shlex.quote(str(bundle / "runtime/env.sh")))
    lines.extend(["export KILNBIO_PYTHON=" + shlex.quote(str(environment / "bin/python")),
                  "export KILNBIO_RSCRIPT=" + shlex.quote(str(environment / "bin/Rscript")),
                  "export KILNBIO_R_LIBRARY=" + shlex.quote(str(environment / "bladder-R")),
                  "export R_LIBS_USER=" + shlex.quote(str(environment / "bladder-R") + os.pathsep + str(r_library)),
                  "exec " + shlex.join([str(x) for x in run] + [str(code_prefix / "bin/kilnbio")]) + ' "$@"'])
    wrapper.write_text("\n".join(lines) + "\n")
    wrapper.chmod(0o755)
    example = smoke(wrapper)
    if source_digest() != source_hash:
        raise ValueError("Source checkout changed during setup; rerun from the final source version")
    bindir = prefix / "bin"
    bindir.mkdir(exist_ok=True)
    pending = bindir / "kilnbio.new"
    pending.write_text('#!/usr/bin/env bash\nexec ' + shlex.quote(str(wrapper)) + ' "$@"\n')
    pending.chmod(0o755)
    pending.replace(bindir / "kilnbio")
    receipt = dict(spec, installed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   bundle=str(bundle), container=args.container, smoke=example)
    pending = prefix / "installation.json.new"
    pending.write_text(json.dumps(receipt, indent=2) + "\n")
    pending.replace(prefix / "installation.json")
    print("Synthetic check passed: 4 donors; 3 with the baseline assay; 2 with both assays.")
    print("Installed. Run " + shlex.quote(str(bindir / "kilnbio")) + " --help")
    print("Resolved package versions: " + str(bundle))
    if "sequencing" in profiles:
        status = container_status(args.container)
        print("Container engine " + args.container + ": " + status)
        if status != "ready":
            print("Sequencing preparation remains: enable the selected container engine. See docs/INSTALLATION.md")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    if sys.version_info < (3, 8):
        parser.error("Python 3.8+ is required to start automatic setup")
    parser.add_argument("--prefix", type=Path, default=Path.home() / ".local/share/kilnbio")
    parser.add_argument("--profile", action="append", choices=[*PROFILES, "all"], help="Repeat to combine profiles; default core")
    parser.add_argument("--container", choices=["docker", "apptainer", "singularity"], default="docker")
    parser.add_argument("--plan", action="store_true", help="Show dependencies without downloading or installing")
    parser.add_argument("--check", action="store_true", help="Recheck an existing installation and its container engine")
    args = parser.parse_args()
    selected = profiles_requested(args.profile)
    try:
        if args.plan:
            print(json.dumps({"profiles": selected, "conda": sorted({p for name in selected for p in PROFILES[name]}),
                              "pip_requirements": [REQUIREMENTS[name] for name in selected if name in REQUIREMENTS],
                              "r_packages": "limma 3.68.5 + statmod 1.5.2; Bioconductor setup_r.R for analysis" if "bladder" in selected else None,
                              "sequencing_runtime": "Java 21 + Nextflow 25.10.4" if "sequencing" in selected else None}, indent=2))
        elif args.check:
            check(args.prefix.expanduser().resolve())
        else:
            install(args, selected)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        parser.exit(1, "KilnBio setup failed: " + str(error) + "\nRerun the same command after fixing the reported issue.\n")


if __name__ == "__main__":
    main()
