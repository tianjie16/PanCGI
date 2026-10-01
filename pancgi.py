import argparse
import json
import subprocess
import sys
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    version = json.loads((root/'RELEASE_MANIFEST.json').read_text())['version']
    commands = {
        'prepare': [sys.executable, str(root/'pancgi_preparation.py'), 'prepare'],
        'run': ['bash', str(root/'scripts/run_pancgi.sh')],
        'inventory-gfa': [sys.executable, str(root/'pancgi_mapping.py'), 'inventory-gfa'],
        'inventory-hal': [sys.executable, str(root/'pancgi_mapping.py'), 'inventory-hal'],
        'mapping-template': [sys.executable, str(root/'pancgi_mapping.py'), 'make-mapping-template'],
        'validate-inputs': [sys.executable, str(root/'pancgi_contract.py')],
        'validate-results': [sys.executable, str(root/'pancgi_validate_results.py')],
    }
    p = argparse.ArgumentParser(prog='pancgi', description='Build and validate graph-native CpG-island catalogues. Each command accepts --help.')
    p.add_argument('--version', action='version', version=f'PanCGI {version}')
    p.add_argument('command', choices=list(commands), nargs='?')
    if len(sys.argv) > 1 and sys.argv[1] in commands:
        command = sys.argv[1]
        arguments = sys.argv[2:]
        if command == 'run':
            arguments = ['--python', sys.executable, *arguments]
        raise SystemExit(subprocess.call([*commands[command], *arguments]))
    args = p.parse_args()
    if args.command is None:
        p.print_help()


if __name__ == '__main__':
    main()
