import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import pancgi_mapping as mapping
from pancgi_contract import digest, file_state


FILES = ('gfa_paths.tsv', 'gfa_paths.tsv.json', 'hal_sequences.tsv', 'hal_genomes.tsv')


def prepare(args):
    output = Path(args.out_dir).resolve()
    if output.exists():
        raise FileExistsError(output)
    states = {str(Path(p).resolve()): file_state(p) for p in (args.gfa, args.hal)}
    output.mkdir(parents=True)
    mapping.inventory_gfa(args.gfa, str(output / 'gfa_paths.tsv'), str(output / 'inventory_work'))
    mapping.inventory_hal(SimpleNamespace(hal=args.hal, output=str(output / 'hal_sequences.tsv'),
        genome_output=str(output / 'hal_genomes.tsv'), docker_bin=args.docker_bin,
        docker_image=args.docker_image, hal_stats=args.hal_stats))
    mapping.make_public_templates(SimpleNamespace(gfa_inventory=str(output / 'gfa_paths.tsv'),
        hal_inventory=str(output / 'hal_sequences.tsv'), out_dir=str(output / 'mapping')))
    if any(file_state(p) != state for p, state in states.items()):
        raise RuntimeError('GFA or HAL changed during preparation')
    report = dict(schema_version=1, inventory_schema_version=2, status='complete',
        sources={name: str(Path(getattr(args, name)).resolve()) for name in ('gfa', 'hal')},
        source_states=states, files={name: digest(output / name) for name in FILES})
    (output / 'preparation.json').write_text(json.dumps(report, indent=2) + '\n')
    return dict(status='complete', output=str(output), mapping_review_required=True)


def validate_and_copy(args):
    source = Path(args.prepared).resolve()
    receipt = (source / 'preparation.json').read_bytes()
    receipt_sha256 = hashlib.sha256(receipt).hexdigest()
    report = json.loads(receipt)
    if report.get('schema_version') != 1 or report.get('inventory_schema_version') != 2 or report.get('status') != 'complete':
        raise ValueError('Unsupported or incomplete preparation')
    expected = {name: str(Path(getattr(args, name)).resolve()) for name in ('gfa', 'hal')}
    if report['sources'] != expected or set(report['source_states']) != set(expected.values()):
        raise ValueError('Preparation does not belong to these GFA/HAL inputs')
    run_states = {p: file_state(p) for p in report['source_states']}
    if any(not isinstance(state, list) or len(state) != 5 or
           any(type(value) is not int for value in state) or
           run_states[p][1:] != state[1:] for p, state in report['source_states'].items()):
        raise ValueError('GFA/HAL changed after preparation; explicitly prepare again')
    if set(report['files']) != set(FILES):
        raise ValueError('Incomplete preparation inventory')
    for name in FILES:
        if digest(source / name) != report['files'][name]:
            raise ValueError(f'Preparation inventory changed: {name}')
    metadata = json.loads((source / 'gfa_paths.tsv.json').read_text())
    if metadata['schema_version'] != 2:
        raise ValueError('Unsupported GFA inventory schema')
    mapping.check_source_fingerprint(args.gfa, metadata['source_fingerprint'], cross_node=True)
    destination = Path(args.out_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    for name in (*FILES, 'preparation.json', 'preparation_reuse.json'):
        if (destination / name).exists():
            raise FileExistsError(destination / name)
    for name in (*FILES, 'preparation.json'):
        shutil.copyfile(source / name, destination / name)
    if (any(digest(destination / name) != report['files'][name] for name in FILES) or
            digest(destination / 'preparation.json') != receipt_sha256):
        raise RuntimeError('Preparation changed while copying')
    if any(file_state(p) != state for p, state in run_states.items()):
        raise RuntimeError('GFA/HAL changed during preparation reuse')
    reuse = dict(status='pass', reused_inventory=True, gfa_source_scans=0, hal_sequence_exports=0,
        source_preparation=str(source), preparation_sha256=receipt_sha256,
        inventory_sha256=report['files'], source_state_fields=['device', 'inode', 'size', 'mtime_ns', 'ctime_ns'],
        cross_node_compared_fields=['path', 'inode', 'size', 'mtime_ns', 'ctime_ns'],
        device_recorded_not_compared_across_nodes=True, within_run_full_state_checked=True,
        original_source_states=report['source_states'], run_source_states=run_states)
    mapping.write_json(destination / 'preparation_reuse.json', reuse)
    return reuse


def main():
    parser = argparse.ArgumentParser(description='Prepare GFA/HAL inventories once, then supply a reviewed mapping to run')
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('prepare', 'validate'):
        p = sub.add_parser(name)
        for option in ('gfa', 'hal', 'out-dir'):
            p.add_argument('--' + option, required=True)
        if name == 'validate':
            p.add_argument('--prepared', required=True)
            p.set_defaults(func=validate_and_copy)
        else:
            p.add_argument('--hal-runtime', choices=['native', 'docker'], default='docker')
            p.add_argument('--docker-bin', default='docker')
            p.add_argument('--docker-image', default='')
            p.add_argument('--hal-stats', default='halStats')
            p.set_defaults(func=prepare)
    args = parser.parse_args()
    if args.command == 'prepare':
        os.environ['PANCGI_HAL_RUNTIME'] = args.hal_runtime
    print(json.dumps(args.func(args), indent=2))


if __name__ == '__main__':
    main()
