"""Run a single instance or a collection with Python in a SageMath environment."""
import argparse
import json
from pathlib import Path
from sage.all import GF, matrix
from minrank_hints import generate, experiment
from progress import activity
from reporting import save_json, print_summary, print_suite_summary


def arguments():
    """Parse the public input, model selection and output paths."""
    parser = argparse.ArgumentParser(description='Support Minors experiments with hints.')
    source = parser.add_mutually_exclusive_group()
    source.add_argument('--config', help='Generator or collection JSON; default: examples/small.json')
    source.add_argument('--instance', help='Saved JSON containing q, r, matrices, hints')
    parser.add_argument('--mode', choices=['both','direct','matching','ordinary','all'], default='both',
                        help='both compares direct substitution and matching')
    parser.add_argument('--max-degree', type=int, help='Largest degree to test')
    parser.add_argument('--no-full-column-reduction', action='store_true')
    parser.add_argument('--save-instance', help='Save a single public instance and hints before solving')
    parser.add_argument('--output', help='JSON path; default: results/results.json or results/small_example2.json')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print each computation phase and its result')
    return parser.parse_args()


def public_input(data):
    """Read saved public matrices, or generate them once from a configuration.

    Saved prime-field instances use integer entries. Hints index Cprime from
    zero. Repeated coordinates are rejected instead of silently overwritten.
    """
    F = GF(data['q'])
    if 'matrices' in data:
        M = [matrix(F, [[F(v) for v in row] for row in A]) for A in data['matrices']]
        rows = data.get('hints', [])
        hints = {(int(i),int(j)): F(value) for i,j,value in rows}
        if len(hints) != len(rows):
            raise ValueError('Duplicate hint coordinates')
        return M, hints
    options = {k:data[k] for k in ('q','m','n','K','r','seed','positions') if k in data}
    return generate(**options)


def save_instance(path, M, r, hints):
    """Save the actual prime-field matrices and hints for later paired comparisons."""
    F = M[0].base_ring()
    if not F.is_prime_field():
        raise ValueError('JSON instance export currently requires a prime field')
    saved = {'q':int(F.order()), 'r':r,
             'matrices':[[[int(a) for a in row] for row in A.rows()] for A in M],
             'hints':[[i,j,int(v)] for (i,j),v in sorted(hints.items())]}
    save_json(path, saved)
    print(f'Instance saved to {path}', flush=True)


def run_one(data, args, output_path, max_degree):
    """Run all selected modes on common inputs and save each completed degree.

    Ctrl+C keeps prior measurements, marks the active run interrupted and exits
    with status 130. A saved partial result cannot resume an unfinished kernel.
    """
    with activity('Loading the public instance' if 'matrices' in data else 'Generating the public instance',
                  args.verbose):
        M, hints = public_input(data)
    if args.save_instance:
        save_instance(args.save_instance, M, data['r'], hints)
    modes = {'both':('direct','matching'), 'all':('ordinary','direct','matching')}.get(args.mode,(args.mode,))
    latest = {}

    def checkpoint(result):
        """Preserve the active phase and complete measurements in the JSON file."""
        latest.clear()
        latest.update(result)
        save_json(output_path, result)

    try:
        result = experiment(M, data['r'], hints, modes=modes,
                            max_degree=max_degree,
                            full_columns=not args.no_full_column_reduction,
                            verbose=args.verbose, checkpoint=checkpoint)
    except (KeyboardInterrupt, Exception) as error:
        if latest:
            for run in latest['runs'].values():
                if run.get('status') == 'running':
                    run['status'] = 'interrupted' if isinstance(error,KeyboardInterrupt) else 'error'
                    if not isinstance(error,KeyboardInterrupt):
                        run['error'] = str(error)
            save_json(output_path, latest)
            print_summary(latest)
            print(f'Partial results saved to {output_path}', flush=True)
        if isinstance(error, KeyboardInterrupt):
            raise SystemExit(130)
        raise
    print_summary(result)
    print(f'Results saved to {output_path}', flush=True)
    return result


def main():
    """Run one experiment, or every explicitly listed instance in a collection.

    Collection paths are relative to their configuration file. Each trial gets
    its own result file as well as an entry in the collection summary. Existing
    input matrices are reused, so direct and matching modes see identical inputs.
    """
    args = arguments()
    if args.max_degree is not None and args.max_degree < 1:
        raise ValueError('max-degree must be positive')
    source = Path(args.instance or args.config or 'examples/small.json')
    data = json.loads(source.read_text())
    if 'instances' not in data:
        output_path = Path(args.output or 'results/results.json')
        degree = args.max_degree if args.max_degree is not None else data.get('max_degree',2)
        run_one(data, args, output_path, degree)
        return
    if args.save_instance:
        raise ValueError('Collection inputs are already saved; --save-instance is for single runs')
    output_path = Path(args.output or 'results/small_example2.json')
    trial_folder = output_path.parent/output_path.stem
    collection = {'description':data.get('description',''), 'trials':[]}
    save_json(output_path, collection)
    total = len(data['instances'])
    try:
        for i, relative in enumerate(data['instances'],1):
            print(f'\nExample {i}/{total}: {relative}', flush=True)
            trial_path = source.parent/relative
            trial = json.loads(trial_path.read_text())
            # Relative subdirectories preserve the grouping by parameter set.
            identity = f'{trial_path.parent.name}-{trial_path.stem}'
            target = trial_folder/(identity+'.json')
            degree = args.max_degree if args.max_degree is not None else data.get('max_degree',3)
            result = run_one(trial, args, target, degree)
            collection['trials'].append({'instance':str(trial_path),
                                         'result_file':str(target), 'result':result})
            save_json(output_path,collection)
    finally:
        save_json(output_path,collection)
    print_suite_summary(collection)
    print(f'Collection results saved to {output_path}', flush=True)


if __name__ == '__main__':
    main()
