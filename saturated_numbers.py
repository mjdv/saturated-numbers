from sympy import isprime

def compute_ds_table():
    """
    Computes the global ds_table, which is used by next to skip
    non-downward-saturated bases. In a function, since it needs to be
    recomputed each time a prime is added to the list primes.
    """
    global ds_table
    ds_table = []
    for p in primes:
        ds_list = []
        acc = p - 1
        for q in primes:
            ct = 0
            while acc % q == 0:
                acc //= q
                ct += 1
            ds_list.append(ct)
        ds_table.append(ds_list)


def divisors_from_fac(fac):
    """
    Given a prime factorization fac (a dictionary of the form prime :
    exponent), returns a generator of all divisors of the factorized 
    number.
    """
    if len(fac) == 0:
        yield 1
        return
    for p in fac: break
    rem_fac = {q : fac[q] for q in fac if q != p}
    for div in divisors_from_fac(rem_fac):
        for e in range(fac[p] + 1):
            yield div*p**e

def num_divisors_from_fac(fac):
    """
    Given a prime factorization fac (a dictionary of the form prime :
    exponent), returns the number of divisors of the factorized number.
    """
    acc = 1
    for p in fac:
        acc *= fac[p] + 1
    return acc

def next(base_value, base_exps):
    """
    Uses global lists primes, mins, ds_table and global cutoff to determine
    the next valid base_value (in lexicographic ordering). That is: a downward
    saturated base that is below the cutoff.

    If there is a next valid base, returns new base_value, base_exps --
    base_exps is also changed in-place.

    If there is no next valid base, returns 0, [].
    """
    for idx in reversed(range(len(primes))):
        if base_exps[idx] == 0 and \
                not all(ce >= re for ce, re in zip(base_exps, ds_table[idx])):
            continue
        base_value *= primes[idx]
        base_exps[idx] += 1
        if base_value < cutoff:
            return base_value, base_exps
        while base_exps[idx] > mins[idx]:
            base_exps[idx] -= 1
            base_value //= primes[idx]
    return 0, []

def next_nondominating(base_value, base_exps):
    """
    Uses global lists primes, mins and global cutoff to determine the next
    valid base in lexicographic ordering which is not a multiple of this base
    (i.e., at least as large on every coordinate). Note that this base is not
    guaranteed to be downward saturated.

    If there is such a base, returns base_value, base_exps -- base_exps is also changed
    in-place.

    If there is no such base, returns 0, [].
    """
    idx = len(primes) - 1
    while idx > 0 and base_exps[idx] == mins[idx]:
        idx -= 1
    if idx == 0:
        return 0, []
    while base_exps[idx] > mins[idx]:
        base_exps[idx] -= 1
        base_value //= primes[idx]
    base_exps[idx-1] += 1
    base_value *= primes[idx-1]
    return base_value, base_exps


def saturated(n, fac):
    """
    Returns whether n with factorization fac is saturated.
    """
    if not downward_saturated(n, fac): return False

    for div in divisors_from_fac(fac):
        if (div + 1) not in fac and isprime(div + 1):
            return False
    return True

def downward_saturated(n, fac):
    """
    Returns whether n with factorization fac is downward saturated.
    """
    for prime in fac:
        if n % (prime - 1) != 0:
            return False
    return True

def dominates_exceeding_base(base_exps):
    """
    Uses global list exceeded_bases to check if base_exps dominates any of them (and
    thus necessarily also will exceed the cutoff while saturating).
    """
    for idx, etup in enumerate(exceeded_bases):
        if all(i >= ei for i, ei in zip(base_exps, etup)):
            break
    else:
        return False
    for j in reversed(range(idx)):
        exceeded_bases[j+1] = exceeded_bases[j]
    exceeded_bases[0] = etup
    return True

# primes contains all currently known primes which divide a saturation number;
# this list is expanded as the program runs.
primes = [2, 3, 7, 43]

# Table used to check whether a given base is downward saturated, and used by
# function next to skip non-ds bases.
ds_table = []
compute_ds_table()

# The program lists all saturated numbers up to cutoff.
cutoff = 10**1000

# Minimal exponent of each prime. Note that each saturated number is divided by
# 2 * 3 * 7 * 43 = 1806.
mins = [1, 1, 1, 1]

# base_value is the current base we are looking to saturate; its prime
# factorization is encoded in base_exps.
base_exps = mins[:]
base_value = 1
for p, e in zip(primes, base_exps):
    base_value *= p**e

# result tracks all saturated numbers (and their factorizations) found.
result = []

# exceeded_bases tracks those bases (as lists of exponents) whose saturation
# sequence exceeds the cutoff; we use it to skip bases which are multiples of
# those.
#
# We seed it with (1, 1, 4, 1) = 2 * 3 * 7^4 * 43. For large cutoffs, it takes
# very long to verify that this number does not divide a saturated number below
# the cutoff. It has been verified (by me) up to 1e1000. If you want to verify
# independently or search a larger cutoff, replace exceeded_bases with the
# empty list.
exceeded_bases = [(1, 1, 4, 1)]

# saturated_gaps tracks those bases n and their saturations n'; then if we have
# a base such that n | m | n', we can skip it.
saturated_gaps = []

while base_value:
    # Check if a divisor of this base already exceeded the cutoff. Then there is
    # no reason to check this base; it must exceed the cutoff. We then move on
    # to the next base which is not a multiple of this base.
    if dominates_exceeding_base(base_exps):
        base_value, base_exps = next_nondominating(base_value, base_exps)
        continue

    fac = {p : base_exps[i] for i, p in enumerate(primes) if base_exps[i]}

    # If this base is not downward saturated, we can skip it; we will encounter
    # the downward saturated base in a later iteration anyway.
    if not downward_saturated(base_value, fac):
        base_value, base_exps = next(base_value, base_exps)
        continue

    # If we have already seen a base n which saturated into n' such that
    # n | base | n', we already know base will saturate into n', and we
    # can skip it.
    in_saturated_gap = False
    for low, high in saturated_gaps:
        if all(ei <= ej <= ek for ei, ej, ek in zip(low, base_exps, high)):
            in_saturated_gap = True
            break

    if in_saturated_gap:
        base_value, base_exps = next(base_value, base_exps)
        continue

    saturating_value = base_value
    while True:
        previous_value = saturating_value

        # We add primes to the base one at a time, to avoid useless searching
        # when we have already passed the cutoff.
        new_prime = 0
        for div in divisors_from_fac(fac):
            if isprime(div + 1) and saturating_value % (div + 1) != 0:
                saturating_value = saturating_value * (div + 1)
                fac[div + 1] = 1
                break

        # We have passed the cutoff; break.
        if saturating_value > cutoff:
            exceeded_bases.append(tuple(base_exps))
            print(f"{base_exps} exceeds, #{len(exceeded_bases)}; {len(fac)} prime divisors, {num_divisors_from_fac(fac)} divisors")
            break

        # We have found a saturating number.
        if saturating_value == previous_value:
            print(f"Found: {saturating_value}, {fac}.")
            for p in fac:
                # If any new primes have been found, add them to the primes
                # list and update the other data accordingly.
                if p not in primes:
                    print(f"NEW PRIME: {p}")
                    primes.append(p)
                    base_exps.append(0)
                    mins.append(0)
                    exceeded_bases = [ex_basis + (0,) for ex_basis in exceeded_bases]
                    saturated_gaps = [(a + (0,), b + (0,)) for a, b in saturated_gaps]
                    compute_ds_table()
            saturated_base = [fac[p] if p in fac else 0 for p in primes]
            saturated_gaps.append((tuple(base_exps), tuple(saturated_base)))
            result.append((saturating_value, fac))
            break

    # Finding the next base to use.
    base_value, base_exps = next(base_value, base_exps)

print("Saturated numbers with their factorizations:")
for value, fac in sorted(result):

    # Uncomment this line to double check at the end that all numbers found are
    # actually saturated. (This is somewhat slow.)
    # assert saturated(value, fac)

    print(value, fac)

print()
print("Saturation primes:")
for p in primes:
    print(p)

print()
print("Saturated numbers as lists of exponents:")
for _, fac in sorted(result):
    print([fac[p] if p in fac else 0 for p in primes])

print()
print(f"{len(result)} saturated numbers found, {len(primes)} saturation primes found.")
