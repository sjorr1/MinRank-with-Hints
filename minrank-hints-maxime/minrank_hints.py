"""Support Minors experiments using SageMath. Indices are zero-based."""
from itertools import combinations
from math import comb
from time import perf_counter
from sage.all import GF, matrix, vector, identity_matrix, zero_matrix, set_random_seed
from progress import activity


def monomials(d, degree):
    """Yield exponent tuples of the requested total degree in d variables."""
    if d == 0:
        if degree == 0:
            yield ()
        return
    if d == 1:
        yield (degree,)
        return
    for a in range(degree + 1):
        for tail in monomials(d - 1, degree - a):
            yield (a,) + tail


def perfect_matching(T, r, mask):
    """Return a row-to-column matching, or None if no perfect matching exists.

    Identity columns have their forced edges. Other entries are allowed
    unless their positions belong to the normalized hint mask.
    """
    neighbours = {i: [j for j in T if (i == j if j < r
                       else (i, j-r) not in mask)] for i in range(r)}
    owner = {}

    def augment(i, seen):
        for j in neighbours[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in owner or augment(owner[j], seen):
                owner[j] = i
                return True
        return False

    for i in range(r):
        if not augment(i, set()):
            return None
    return {i: j for j, i in owner.items()}


def elementary_set(r, i, j):
    """Return the minor replacing identity column i by column r+j; indices start at zero."""
    return tuple(k for k in range(r) if k != i) + (r+j,)


def generate(q=101, m=7, n=10, K=9, r=5, seed=1, positions=None):
    """Generate independent public matrices and hints for a planted instance.

    Work over GF(q). The K public matrices have size m by n; a random coefficient
    vector with first entry one gives H*C, with rank(H)=r and C=(I_r | Cprime).
    Hint positions are zero-based Cprime indices. Save the actual returned
    matrices to compare methods on identical inputs, independently of seeds.
    """
    F = GF(q)
    if not (1 <= r < min(m, n) and 1 <= K <= m*n):
        raise ValueError('Require 1 <= r < min(m,n) and 1 <= K <= m*n')
    set_random_seed(seed)
    C = identity_matrix(F, r).augment(matrix(F, r, n-r,
                                  [F.random_element() for _ in range(r*(n-r))]))
    while True:
        H = matrix(F, m, r, [F.random_element() for _ in range(m*r)])
        if H.rank() == r:
            break
    E = H*C
    while True:
        x = vector(F, [F.one()] + [F.random_element() for _ in range(K-1)])
        M = [zero_matrix(F, m, n)] + [matrix(F, m, n,
                 [F.random_element() for _ in range(m*n)]) for _ in range(K-1)]
        M[0] = E - sum((x[k]*M[k] for k in range(1, K)), zero_matrix(F, m, n))
        if matrix(F, [A.list() for A in M]).rank() == K:
            break
    if positions is None:
        positions = [(i, j) for i in range(min(3, r))
                     for j in range(max(0, n-r-3), n-r)]
    hints = {(i, j): C[i, r+j] for i, j in positions}
    return M, hints


def reduce_full_columns(M, r, hints):
    """Eliminate the actual constraints imposed by fully revealed Cprime columns.

    Return reduced matrices, remaining hints with new column indices, a basis B,
    the deleted Cprime columns, and the constraint rank. A reduced coefficient
    row vector u gives the original vector x=u*B. Independence is not assumed."""
    F = M[0].base_ring()
    m, n, K = M[0].nrows(), M[0].ncols(), len(M)
    full = [j for j in range(n-r) if all((i, j) in hints for i in range(r))]
    constraints = [[A[k, r+j] - sum((hints[i, j]*A[k, i]
                      for i in range(r)), F.zero()) for A in M]
                   for j in full for k in range(m)]
    L = matrix(F, constraints) if constraints else matrix(F, 0, K)
    B = L.right_kernel().basis_matrix()
    keep = list(range(r)) + [r+j for j in range(n-r) if j not in full]
    reduced = [sum((v[k]*M[k] for k in range(K)), zero_matrix(F, m, n))
               .matrix_from_columns(keep) for v in B.rows()]
    rem = [j for j in range(n-r) if j not in full]
    new_hints = {(i, rem.index(j)): value for (i, j), value in hints.items()
                 if j in rem}
    return reduced, new_hints, B, full, int(L.rank())


def model(M, r, hints, mode):
    """Return coordinates and sparse bilinear equations.

    direct: substitute the actual hinted elementary minors.
    matching: first normalize the hints, then remove all forbidden minors.
    ordinary: retain every minor and ignore the remaining hints.
    """
    F = M[0].base_ring()
    m, n, d = M[0].nrows(), M[0].ncols(), len(M)
    all_T = list(combinations(range(n), r))
    T0 = tuple(range(r))
    removed = {}
    if mode == 'matching':
        P = identity_matrix(F, n)
        for (i, j), value in hints.items():
            P[i, r+j] = -value
        M = [A*P for A in M]
        Ts = [T for T in all_T if perfect_matching(T, r, set(hints)) is not None]
    elif mode == 'direct':
        removed = {elementary_set(r, i, j): (-1)**(r-1-i)*value
                   for (i, j), value in hints.items()}
        Ts = [T for T in all_T if T not in removed]
    elif mode == 'ordinary':
        Ts = all_T
    else:
        raise ValueError('mode must be ordinary, direct, or matching')
    index = {T: t for t, T in enumerate(Ts)}
    equations, blocks = [], 0
    for L in combinations(range(n), r+1):
        terms = []
        for a, j in enumerate(L):
            T = L[:a] + L[a+1:]
            if T in index:
                terms.append((j, index[T], F((-1)**a)))
            elif T in removed:
                terms.append((j, index[T0], F((-1)**a)*removed[T]))
        if not terms:
            continue
        blocks += 1
        for u in range(m):
            equation = {}
            for j, t, sign in terms:
                for k in range(d):
                    key = (k, t)
                    equation[key] = equation.get(key, F.zero()) + sign*M[k][u, j]
            equations.append({key: a for key, a in equation.items() if a})
    return Ts, equations, blocks


def macaulay(F, d, s, equations, b):
    """Build the degree-(b,1) matrix and its exponent-to-column-block map.

    Multiply each bilinear equation by every degree-(b-1) x-monomial. Columns
    are ordered by exponent tuple, then minor index. Preserve zero rows and
    impose neither field equations nor additional determinantal identities.
    """
    upper, lower = list(monomials(d, b)), list(monomials(d, b-1))
    rows, cols = len(equations)*len(lower), s*len(upper)
    position = {v: a for a, v in enumerate(upper)}
    entries = {}
    row = 0
    for v in lower:
        for equation in equations:
            for (k, t), value in equation.items():
                w = list(v)
                w[k] += 1
                col = position[tuple(w)]*s+t
                entries[row, col] = entries.get((row, col), F.zero()) + value
            row += 1
    return matrix(F, rows, cols, entries, sparse=True), position


def validate(M, r, hints, x):
    """Check a nonzero candidate against the original matrices and hint values.

    Require exact rank r and independent first r columns. Recover C from S*C=E
    and compare each revealed Cprime entry with its original hint value."""
    if not any(x):
        return False
    F = M[0].base_ring()
    E = sum((x[k]*A for k, A in enumerate(M)), zero_matrix(F, M[0].nrows(), M[0].ncols()))
    S = E.matrix_from_columns(range(r))
    if E.rank() != r or S.rank() != r:
        return False
    C = S.solve_right(E)
    return S*C == E and all(C[i, r+j] == value for (i, j), value in hints.items())


def experiment(M, r, hints, modes=('direct', 'matching'), max_degree=2,
               full_columns=True, verbose=False,
               checkpoint=None):
    """Compare measured ranks and validated recovery on one common instance.

    M contains SageMath matrices over the same finite field; hints maps zero-based
    Cprime positions to field elements. Ordinary mode ignores all hints. The
    other modes optionally eliminate fully known columns using their actual rank.
    A one-dimensional kernel is decoded and validated in the original instance.
    A larger kernel is not evidence that no solution exists.

    verbose prints each phase before it starts and after it finishes.
    checkpoint(result), if provided, is
    called at phase boundaries and after each completed degree to preserve work.
    No field equations are imposed. Return a JSON-ready dictionary.
    """
    if max_degree < 1:
        raise ValueError('max_degree must be positive')
    if not M:
        raise ValueError('Empty public matrix list')
    F = M[0].base_ring()
    m, n, K = M[0].nrows(), M[0].ncols(), len(M)
    if not 1 <= r <= min(m, n) or any(A.base_ring() != F or A.dimensions() != (m,n) for A in M):
        raise ValueError('Inconsistent matrix dimensions or fields')
    if any(not (0 <= i < r and 0 <= j < n-r) for i,j in hints):
        raise ValueError('Hint coordinates must index C prime, starting at zero')
    hints = {key: F(value) for key,value in hints.items()}
    output = {'q': int(F.order()), 'm': m, 'n': n, 'K': K, 'r': r, 'runs': {}}
    def save():
        if checkpoint is not None:
            checkpoint(output)

    for mode in modes:
        run = {'status': 'running', 'phase': 'preprocessing', 'attempts': []}
        output['runs'][mode] = run
        save()
        if verbose:
            print(f'\n[{mode}] Preprocessing the instance...', flush=True)
        start = perf_counter()
        with activity(f'[{mode}] Reducing fully known columns', verbose):
            # Ordinary means no hint information, including no full-column reduction.
            if full_columns and mode != 'ordinary':
                A, H, B, full, rho = reduce_full_columns(M, r, hints)
            else:
                A, H, B, full, rho = M, hints, identity_matrix(F, K), [], 0
        d = len(A)
        run.update(d=d, n0=n-len(full), full_columns=full, linear_constraint_rank=rho)
        if not d:
            run['status'] = 'no_nonzero_coefficient_vector'
            save()
            continue
        with activity(f'[{mode}] Constructing equations', verbose):
            Ts, equations, blocks = model(A, r, H, mode)
        s, T0 = len(Ts), tuple(range(r))
        run.update(s=s, blocks=blocks, preprocessing_seconds=perf_counter()-start)
        if verbose:
            print(f'[{mode}] d={d}; n0={n-len(full)}; {s} minor coordinates; '
                  f'{blocks} equation blocks.', flush=True)
        save()
        for b in range(1, max_degree+1):
            rows = len(equations)*comb(d+b-2,b-1)
            cols = s*comb(d+b-1,b)
            run.update(phase='matrix construction', current_degree=b)
            save()
            tick = perf_counter()
            with activity(f'[{mode}] Degree {b}: building {rows:,} x {cols:,} matrix',
                          verbose):
                Mb, position = macaulay(F, d, s, equations, b)
            built = perf_counter()
            run['phase'] = 'right kernel'
            save()
            with activity(f'[{mode}] Degree {b}: computing the right kernel',
                          verbose):
                kernel = Mb.right_kernel()
            solved = perf_counter()
            record = {'b': b, 'rows': Mb.nrows(), 'columns': Mb.ncols(),
                      'rank': int(Mb.ncols()-kernel.dimension()),
                      'kernel_dimension': int(kernel.dimension()),
                      'build_seconds': built-tick, 'solve_seconds': solved-built}
            run['attempts'].append(record)
            run['phase'] = 'validation' if kernel.dimension() == 1 else 'degree complete'
            if verbose:
                print(f'[{mode}] Degree {b}: rank {record["rank"]:,}; '
                      f'kernel dimension {record["kernel_dimension"]:,}.', flush=True)
            save()
            if kernel.dimension() == 1:
                v = kernel.basis()[0]
                t0 = Ts.index(T0)
                for k in range(d):
                    w = [0]*d
                    w[k] = b
                    z = v[position[tuple(w)]*s+t0]
                    if not z:
                        continue
                    ratios = []
                    for j in range(d):
                        w = [0]*d
                        w[k] = b-1
                        w[j] += 1
                        ratios.append(v[position[tuple(w)]*s+t0]/z)
                    x = vector(F, ratios)*B
                    with activity(f'[{mode}] Degree {b}: validating the candidate',
                                  verbose):
                        good = validate(M, r, hints if mode != 'ordinary' else {}, x)
                    record['valid_solution'] = bool(good)
                    record['x'] = [str(a) for a in x]
                    break
                if record.get('valid_solution'):
                    run['status'] = 'recovered'
                    break
            save()
        if run['status'] == 'running':
            run['status'] = 'not_recovered_in_tested_range'
        run['phase'] = 'finished'
        run['total_seconds'] = perf_counter()-start
        save()
    return output


def ordinary_count(b, m, n, r, d):
    """The ordinary alternating count; not a certified hinted rank."""
    choose = lambda a, k: comb(a, k) if 0 <= k <= a else 0
    return sum((-1)**(i+1)*choose(n,r+i)*choose(m+i-1,i)*choose(d+b-i-1,b-i)
               for i in range(1,b+1))
