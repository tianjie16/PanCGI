import argparse
import csv
import gzip
import json
import multiprocessing
import os
import tempfile
import time
import unittest
from unittest.mock import patch
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import cpgi_nr_prod as prod
import pancgi_feature_build as feature_build
import pancgi_features as feature_store
import pancgi_graph_unfold as unfold
import pancgi_mapping as mapping
from pancgi_parallel import completed_tasks


def write_tsv(path, columns, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def payload(path):
    return gzip.decompress(path.read_bytes()) if str(path).endswith('.gz') else path.read_bytes()


def delayed_task(value):
    time.sleep(0.01 * (5 - value % 5))
    return value, os.getpid()


def killed_task(value):
    if value == 1:
        os._exit(23)
    time.sleep(0.05)
    return value


class ParallelSchedulingTests(unittest.TestCase):
    def test_bounded_submission_and_unique_results(self):
        for workers in (2, 8):
            seen = []
            def tasks():
                for i in range(20):
                    seen.append(i)
                    yield i
            results = {}
            for ordinal, result in completed_tasks(delayed_task, tasks(), workers):
                results[ordinal] = result
                self.assertLessEqual(len(seen), workers + len(results) - 1)
            self.assertEqual([results[i][0] for i in range(20)], list(range(20)))
            self.assertGreater(len({r[1] for r in results.values()}), 1)

    def test_worker_death_propagates_and_joins_children(self):
        with self.assertRaises(BrokenProcessPool):
            list(completed_tasks(killed_task, range(10), 2))
        self.assertEqual(multiprocessing.active_children(), [])


class ParallelGraphTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self, opaque=False, missing=False, bad_link=False):
        a, b = ('opaque.alpha', 'unitig-B') if opaque else ('1', '2')
        text = f'S\t{a}\tAAA\nS\t{b}\tCC\nL\t{b}\t-\t{a}\t-\t' + ('1M\n' if bad_link else '0M\n')
        for i in range(12):
            text += f'W\tindividual-{i}\tcustom-hap\tcontig-{i}\t0\t5\t>{a}>{b}\n'
            text += f'W\tindividual-{i}\tcustom-hap\tcontig-{i}\t5\t8\t<' + ('9999' if missing and i == 4 else a) + '\n'
        text += f'P\tpath name without delimiters\t{a}+,{b}+,{a}-\t0M,0M\n'
        text += f'P\tlinked-path\t{a}+,{b}+\t*\n'
        gfa, inv, contigs = (self.root / name for name in ('graph.gfa', 'inventory.tsv', 'contigs.tsv'))
        gfa.write_text(text)
        mapping.inventory_gfa(str(gfa), str(inv), str(self.root / 'inventory'))
        rows = mapping.read_tsv(str(inv), mapping.GFA_INVENTORY_COLUMNS)
        mapping.write_tsv(contigs, [dict(genome_id=f'g{i}', contig_id=f'c{i}', gfa_path_id=row['gfa_path_id'],
            mapping_status='confirmed') for i, row in enumerate(rows)], mapping.CONTIG_COLUMNS)
        return gfa, inv, contigs, {row['gfa_path_id']: 5 if row['raw_path_name'] == 'linked-path' else 8 for row in rows}

    def run_graph(self, inputs, workers, name, compression='gzip'):
        gfa, inv, contigs, lengths = inputs
        out = self.root / name
        result = unfold.unfold_gfa(str(gfa), str(out), contigs=str(contigs), gfa_inventory=str(inv),
            expected_lengths=lengths, threads=workers, compression=compression)
        return out, result

    def test_dense_and_opaque_exact_at_one_two_eight_workers(self):
        for opaque in (False, True):
            with self.subTest(opaque=opaque):
                case = self.root / ('opaque' if opaque else 'dense')
                case.mkdir()
                original_root, self.root = self.root, case
                inputs = self.prepare(opaque=opaque)
                first, serial = self.run_graph(inputs, 1, 'one')
                for workers in (2, 8):
                    other, result = self.run_graph(inputs, workers, f'w{workers}')
                    self.assertEqual((first / 'pathbed_inventory.tsv').read_bytes(), (other / 'pathbed_inventory.tsv').read_bytes())
                    for path in first.glob('*.bed.gz'):
                        self.assertEqual(path.read_bytes(), (other / path.name).read_bytes())
                    for key in ('n_pathbed_rows', 'pathbed_span_bp', 'n_confirmed_paths', 'source_read_passes', 'segment_index_builds', 'link_index_builds'):
                        self.assertEqual(serial[key], result[key])
                    self.assertEqual(result['workers_used'], workers)
                    events = [json.loads(x) for x in (other / 'progress.jsonl').read_text().splitlines()]
                    passed = [x['gfa_path_id'] for x in events if x['phase'] == 'validation' and x['event'] == 'path_complete']
                    self.assertEqual(len(passed), len(set(passed)))
                    self.assertEqual(set(passed), set(inputs[3]))
                self.root = original_root

    def test_failure_does_not_publish_paths_or_summary(self):
        inputs = self.prepare(missing=True)
        for workers in (1, 2, 8):
            with self.subTest(workers=workers), self.assertRaisesRegex(ValueError, 'undefined segment'):
                self.run_graph(inputs, workers, f'bad{workers}')
            out = self.root / f'bad{workers}'
            self.assertFalse(list(out.glob('*.bed.gz')))
            self.assertFalse((out / 'pathbed_inventory.tsv').exists())
            self.assertFalse((out / 'unfold_summary.json').exists())

    def test_readonly_link_checks_are_not_skipped(self):
        inputs = self.prepare(opaque=True, bad_link=True)
        with self.assertRaisesRegex(ValueError, 'link_overlap_unsupported'):
            self.run_graph(inputs, 2, 'bad_link')

    def test_expected_length_and_invalid_workers_fail(self):
        inputs = self.prepare()
        inputs[3][next(iter(inputs[3]))] += 1
        with self.assertRaisesRegex(ValueError, 'HAL span mismatch'):
            self.run_graph(inputs, 2, 'wrong_length')
        with self.assertRaisesRegex(ValueError, 'positive integer'):
            self.run_graph(inputs, 0, 'zero')
        self.assertFalse((self.root / 'zero').exists())

    def test_plain_output_and_parent_source_check(self):
        inputs = self.prepare(opaque=True)
        first, _ = self.run_graph(inputs, 1, 'plain1', compression='none')
        other, _ = self.run_graph(inputs, 2, 'plain2', compression='none')
        for path in first.glob('*.bed'):
            self.assertEqual(path.read_bytes(), (other / path.name).read_bytes())
        original = mapping.check_source_fingerprint
        calls = []
        def checked(*args, **kwargs):
            calls.append(1)
            if len(calls) == 4:
                raise ValueError('source changed before publication')
            return original(*args, **kwargs)
        with patch.object(mapping, 'check_source_fingerprint', side_effect=checked):
            with self.assertRaisesRegex(ValueError, 'source changed before publication'):
                self.run_graph(inputs, 2, 'changed')
        self.assertEqual(len(calls), 4)
        self.assertFalse(list((self.root / 'changed').glob('*.bed.gz')))
        self.assertFalse((self.root / 'changed/unfold_summary.json').exists())


class ParallelFeatureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.psl, self.paths = self.root / 'psl', self.root / 'paths'
        self.psl.mkdir()
        self.paths.mkdir()
        self.catalog = self.root / 'catalog.tsv'
        self.rows = []
        inventory = []
        for i in range(12):
            label, contig = f'g{i}', f'c{i}'
            kind = 'reference_primary' if i == 0 else ('reference_comparison' if i == 1 else 'assembly')
            bed, fasta, sv = (self.root / f'{label}.{ext}' for ext in ('bed', 'fa', 'sv.tsv'))
            intervals = [] if i == 11 else [(10, 30), (60, 90), (120, 140)]
            bed.write_text(''.join(f'{label}#0#{contig}\t{s}\t{e}\tCGI\t{e-s}\t2\t12\t20\t60\t1\n' for s, e in intervals))
            fasta.write_text(''.join(f'>{label}#0#{contig}:{s}-{e}\n' + 'CG' * ((e-s)//2) + '\n' for s, e in intervals))
            sv.write_text('#CHROM\tPOS\tVCF_END\tID\tSVTYPE\tcontig\tstart\tend\tSVLEN\n')
            path = self.paths / f'{label}.bed'
            path.write_text(f'{contig}\t0\t40\t>1\n{contig}\t40\t80\t<2\n')
            inventory.append(dict(genome_id=label, contig_id=contig, output_file=path.name))
            psl_rows = []
            for s, e in intervals:
                psl_rows.append([f'{label}#0#{contig}:{s}-{e}', e-s, 0, 0, 0, 0, 0, 0, 0,
                    '++', contig, 200, s, e, 'primary-contig', 1000, s, e, 1, f'{e-s},', f'{s},', f'{s},'])
            with (self.psl / f'{label}.to_primary.psl').open('w') as handle:
                csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(psl_rows)
            self.rows.append(dict(label=label, kind=kind, graph_sample=label, graph_hap='0', hal_genome=label,
                bed=str(bed), cpgi_fa=str(fasta), sv_tsv=str(sv), path_bed_dir=str(self.paths)))
        write_tsv(self.catalog, prod.INTERNAL_CATALOG_COLUMNS, self.rows)
        write_tsv(self.paths / 'pathbed_inventory.tsv', ['genome_id', 'contig_id', 'output_file'], inventory)

    def tearDown(self):
        self.temp.cleanup()

    def arguments(self, workers, name, compressed=True):
        suffix = '.gz' if compressed else ''
        return argparse.Namespace(catalog=str(self.catalog), out=str(self.root / (name + '.jsonl' + suffix)),
            excluded=str(self.root / (name + '.excluded.tsv' + suffix)), hal_psl_dir=str(self.psl),
            min_graph_cov=.95, flank_bp=1000, flank_max_steps=32, hal_min_coverage=.5,
            hal_min_identity=0., sv_contig_coordinate_base=0, threads=workers)

    def test_one_two_eight_workers_preserve_order_exclusions_and_global_weights(self):
        first = self.arguments(1, 'one')
        feature_build.run(first)
        records = [json.loads(x) for x in payload(Path(first.out)).splitlines()]
        self.assertEqual([r['label'] for r in records], [f'g{i}' for i in range(11)])
        self.assertEqual(len(payload(Path(first.excluded)).splitlines()), 23)
        for workers in (2, 8):
            args = self.arguments(workers, f'w{workers}')
            feature_build.run(args)
            self.assertEqual(payload(Path(first.out)), payload(Path(args.out)))
            self.assertEqual(payload(Path(first.excluded)), payload(Path(args.excluded)))
            db = self.root / f'w{workers}.sqlite'
            feature_store.build(args.out, db, 3, 8)
            self.assertEqual(prod.collect_shingle_df(first.out, 3, 8), prod.collect_shingle_df(str(db), 3, 8))
            events = [json.loads(x) for x in Path(args.out + '.progress.jsonl').read_text().splitlines()]
            self.assertEqual(events[-1]['event'], 'complete')
            self.assertEqual(events[-1]['input_records'], 33)
            self.assertEqual(events[-1]['features_written'], 11)
            self.assertEqual(events[-1]['excluded'], 22)
            self.assertEqual([x['ordinal'] for x in events if x['event'] == 'genome_merged'], list(range(12)))

    def test_plain_text_format_preserved(self):
        a, b = self.arguments(1, 'plain1', False), self.arguments(2, 'plain2', False)
        feature_build.run(a)
        feature_build.run(b)
        self.assertEqual(payload(Path(a.out)), payload(Path(b.out)))
        self.assertEqual(payload(Path(a.excluded)), payload(Path(b.excluded)))

    def test_genome_failure_does_not_publish_outputs(self):
        Path(self.rows[4]['cpgi_fa']).write_text('')
        args = self.arguments(2, 'bad')
        with self.assertRaisesRegex(RuntimeError, 'g4.*lacks'):
            feature_build.run(args)
        self.assertFalse(Path(args.out).exists())
        self.assertFalse(Path(args.excluded).exists())
        events = [json.loads(x) for x in Path(args.out + '.progress.jsonl').read_text().splitlines()]
        self.assertEqual(events[-1]['event'], 'error')
        self.assertNotIn('complete', [e['event'] for e in events])

    def test_repeated_run_and_invalid_workers_rejected(self):
        args = self.arguments(1, 'once')
        feature_build.run(args)
        with self.assertRaises(FileExistsError):
            feature_build.run(args)
        with self.assertRaisesRegex(ValueError, 'positive integer'):
            feature_build.run(self.arguments(0, 'zero'))


if __name__ == '__main__':
    unittest.main()
