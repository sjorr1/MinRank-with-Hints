"""Readable terminal summaries and atomic JSON checkpoints; standard Python only."""
import json
from pathlib import Path
import os
import tempfile


def save_json(path, result):
    """Replace the result file atomically, preserving its last complete checkpoint."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=target.name+'.', dir=target.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def print_summary(result):
    """Display dimensions, measured ranks, timings and validation without raw dicts."""
    print('\nMinRank with Hints — results', flush=True)
    print(f"Field: GF({result['q']}); m={result['m']}, n={result['n']}, "
          f"K={result['K']}, r={result['r']}")
    names = {'ordinary': 'Ordinary Support Minors',
             'direct': 'Direct hint substitution', 'matching': 'Matching preprocessing'}
    for mode, run in result['runs'].items():
        print('\n'+names.get(mode, mode))
        print('Status: '+run.get('status', 'unknown').replace('_', ' '))
        if 's' in run:
            print(f"Reduced parameters: d={run['d']}, n0={run['n0']}; "
                  f"{run['s']:,} homogeneous minor coordinates; "
                  f"{run['blocks']:,} equation blocks.")
        header = ('Degree', 'Rows', 'Columns', 'Rank', 'Kernel dim.', 'Build (s)', 'Kernel (s)')
        rows = []
        for record in run['attempts']:
            rows.append((str(record['b']), f"{record['rows']:,}",
                         f"{record['columns']:,}", f"{record['rank']:,}",
                         f"{record['kernel_dimension']:,}",
                         f"{record['build_seconds']:.3f}", f"{record['solve_seconds']:.3f}"))
        if rows:
            widths = [max(len(row[i]) for row in [header]+rows) for i in range(len(header))]
            for row in [header]+rows:
                print('  '.join(value.rjust(width) for value,width in zip(row,widths)))
        else:
            print('No degree computation completed.')
        for record in run['attempts']:
            if 'valid_solution' in record:
                print(f"Degree {record['b']} candidate validation: "
                      + ('passed' if record['valid_solution'] else 'failed'))
                print('x = ('+', '.join(record['x'])+')')
        if run.get('status') == 'recovered':
            degree = next(a['b'] for a in run['attempts'] if a.get('valid_solution'))
            print(f'First validated corank-one recovery: degree {degree}.')
        elif run.get('status') == 'not_recovered_in_tested_range':
            print('No validated corank-one recovery in the tested range; '
                  'this does not establish unsolvability.')
    print(flush=True)


def print_suite_summary(collection):
    """Summarize successful degrees across masks, grouped by parameters and method.

    Only validated recoveries contribute to degree extrema. Trial timings and
    matrix measurements remain in the individual JSON files; these summaries
    are observations over the supplied inputs, not general solving bounds.
    """
    groups = {}
    for trial in collection['trials']:
        result = trial['result']
        params = tuple(result[k] for k in ('q','m','n','K','r'))
        for mode, run in result['runs'].items():
            item = groups.setdefault((params,mode), {'count':0,'degrees':[]})
            item['count'] += 1
            degrees = [a['b'] for a in run['attempts'] if a.get('valid_solution')]
            if degrees:
                item['degrees'].append(min(degrees))
    print('\nSmall example 2 — collection summary')
    header = ('(q,m,n,K,r)', 'Method', 'Recovered', 'Min degree', 'Max degree')
    rows=[]
    for (params,mode),item in groups.items():
        degrees=item['degrees']
        rows.append((str(params),mode,f"{len(degrees)}/{item['count']}",
                     str(min(degrees)) if degrees else '-',str(max(degrees)) if degrees else '-'))
    if not rows:
        print('No trial completed.')
        return
    widths=[max(len(row[i]) for row in [header]+rows) for i in range(len(header))]
    for row in [header]+rows:
        print('  '.join(value.ljust(width) for value,width in zip(row,widths)))
