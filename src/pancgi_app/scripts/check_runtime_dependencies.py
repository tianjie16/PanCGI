from __future__ import annotations

import importlib
import json
import sys


REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "pyarrow": "pyarrow",
    "parasail": "parasail",
    "Bio": "biopython",
    "pywfa": "pywfa",
}


def main() -> None:
    versions = {}
    missing = []
    for module_name, package_name in REQUIRED_MODULES.items():
        try:
            module = importlib.import_module(module_name)
        except Exception as exc:
            missing.append({"module": module_name, "package": package_name, "error": repr(exc)})
            continue
        versions[module_name] = str(getattr(module, "__version__", "available"))

    try:
        from pywfa import WavefrontAligner
    except Exception as exc:
        missing.append({"module": "pywfa.WavefrontAligner", "package": "pywfa", "error": repr(exc)})

    report = {
        "status": "fail" if missing else "pass",
        "python": sys.executable,
        "python_version": sys.version.split()[0],
        "modules": versions,
        "missing": missing,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if missing:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
