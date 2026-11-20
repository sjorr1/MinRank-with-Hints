# MinRank with Hints: SageMath experiments

Here is the implementation of our SupportMinors modeling with hints
from the companion paper. It compares ordinary SupportMinors, the
modeling with direct hint substitution, and the addition of matching
preprocessing.

AI disclosure: ChatGPT-6.1 SOL was used to help generate the code for
the matching algorithm, and to format the examples in JSON so that
they could be tested on different platforms.

To run the experiments, use Python from an environment where SageMath
is installed. You can then simply execute

```bash
python run_experiment.py --config examples/small.json --mode all --verbose
```

It stores the outputs in `results/` in JSON format so that they could be
easily exploited if needed.

## Small example 1

`examples/small.json` has (q,m,n,K,r)=(7,4,5,5,2), with two hints at Cprime
positions (0,0) and (1,1). It generates one planted instance using seed 1.

## Small example 2: the five sets from the second archive

`examples/small_example2.json` lists the 50 inputs corresponding to
Table 4 of the companion paper: five parameter sets, each with ten
different hint masks on the same public instance.

| Set | (q,m,n,K,r) | Number of hints | Trials |
|---|---|---:|---:|
| 1 | (7,4,5,5,2) | 2 | 10 |
| 2 | (7,4,5,6,2) | 3 | 10 |
| 3 | (7,5,6,7,2) | 3 | 10 |
| 4 | (7,5,6,8,2) | 4 | 10 |
| 5 | (7,5,6,7,3) | 5 | 10 |

The public matrices and hints are provided explicitly as JSON files.

Run the same matching experiments:

```bash
python run_experiment.py --config examples/small_example2.json \
  --mode matching --max-degree 3 --verbose \
  --output results/small_example2.json
```

Compare direct substitution and matching on every input:

```bash
python run_experiment.py --config examples/small_example2.json \
  --mode both --max-degree 3 --verbose \
  --output results/small_example2-comparison.json
```

Run just the fifth parameter set:

```bash
python run_experiment.py --config examples/small_example2/set-5.json \
  --mode both --max-degree 3 --verbose \
  --output results/small_example2-set5.json
```

Run an individual input:

```bash
python run_experiment.py \
  --instance examples/small_example2/set-5/trial-01.json \
  --mode both --max-degree 3 --verbose \
  --output results/set5-trial01.json
```

## Larger example:

```bash
python run_experiment.py \
  --config examples/paper-parameters.json --mode matching \
  --max-degree 2 --verbose --output results/matching.json \
  --save-instance results/instance.json

python run_experiment.py \
  --instance results/instance.json --mode direct \
  --max-degree 2 --verbose --output results/direct.json
```

This uses (q,m,n,K,r)=(101,7,10,9,5) with the nine hints in rows 1,2,3 of
global support columns 8,9,10. It generates a planted instance with these
parameters.

| Method | Degree | Rows | Columns |
|---|---:|---:|---:|
| Direct | 1 | 1470 | 2187 |
| Direct | 2 | 13230 | 10935 |
| Matching | 1 | 1085 | 1134 |
| Matching | 2 | 9765 | 5670 |


## Providing Custom Inputs

Hint positions are zero-based indices of Cprime, not of the full C=(I_r|Cprime).
For example, position (i,j) means C[i,r+j]. 

For a saved prime-field input, supply JSON keys `q`, `r`, `matrices` (a list of
matrices) and `hints` (triples [i,j,value]). `--save-instance` exports a single
input.


## Miscellaneous

Use `--verbose` or `-v` for more information while running the
experiments. The JSON file is updated during the experiments to preserve
the results of completed degrees. If the computation is interrupted,
these results remain available.


## Files

- `run_experiment.py`: command-line interface
- `minrank_hints.py`: generation, preprocessing, matrix construction and recovery.
- `progress.py`: phase messages.
- `reporting.py`: terminal tables and JSON checkpoints.
- `examples/`: configurations and saved inputs.
