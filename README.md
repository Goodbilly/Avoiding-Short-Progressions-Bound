# Avoiding Short Progressions in Euclidean Ramsey Theory

**Research code and exact computations for improved Euclidean Ramsey bounds.**

This repository accompanies our work on two-colourings of Euclidean spaces avoiding short monochromatic arithmetic progressions. It provides exact arithmetic checks for the uniform 22-term bound, a 17-term result for irrational squared scale, and an investigation of the 21-term bound.

## The problem and the results

Write

$$
\mathbb E^n\not\longrightarrow (\ell_3,\alpha\ell_M)
$$

if there exists a red–blue colouring of $\mathbb R^n$ with **no red three-term collinear arithmetic progression of consecutive spacing $1$** and **no blue $M$-term collinear arithmetic progression of consecutive spacing $\alpha$**. Here $\alpha>0$ and $n\ge1$.

Our main results and investigations are:

| Number of blue terms | Scope | Status |
| --- | --- | --- |
| $M=22$ | Every $\alpha>0$ and every $n\ge1$ | **Proved** |
| $M=17$ | In particular, every $\alpha$ with irrational $\alpha^2$ | **Proved in the stated special regime** |
| $M=21$ | Every $\alpha>0$ except possibly $\alpha^2\in\{1/4,1/9,2/15,3/4\}$ | **Partial result** |

**The uniform $M=21$ bound is not claimed.** The four displayed squared scales are unresolved by the methods implemented here; they are not known counterexamples.

## Quick start

Requirements: Python 3 (standard library only), a C++17 compiler such as GCC or Clang, and a terminal. The C++ code uses `__int128` and some GNU headers, so GCC on Linux/WSL is a convenient choice. Run all commands from the repository root.

```bash
git clone https://github.com/Goodbilly/Avoiding-Short-Progressions-Bound.git
cd Avoiding-Short-Progressions-Bound
```

### 1. Verify the uniform 22-term result and the 17-term computations

```bash
g++ -O2 -std=c++17 verify_uniform22_with_rational17.cpp -o verify_uniform22
./verify_uniform22
```

This is the main independent verification program. It uses integer/rational arithmetic to check:

- the individual $G_{p,22}(u)$ values used in the proof;
- the inequality $G_{33,22}(u)\le8$ throughout $[9/25,16/25]$, via midpoint subdivision and a Lipschitz bound;
- the remaining **187 reduced rational cases**, divided by the code into those handled by the large interval and those requiring additional coefficients;
- the stated 20-term gap computation and the interval used in the 17-term argument.

The code prints its rational-case counts and reports successful checks. A successful execution supports the *computational steps* of the paper; the mathematical reductions and colouring arguments are explained in the manuscript.

### 2. Verify four explicit 21-term colourings

```bash
g++ -O2 -std=c++17 verify_four.cpp -o verify_four
./verify_four
```

The program checks the red-spacing condition and computes exact circular gaps for:

| $\alpha^2$ | $p$ | $s$ | $R$ | $G_{p,21}(\alpha^2R)$ |
| --- | ---: | ---: | ---: | ---: |
| $2/23$ | 5 | 1 | $1$ | $203/207$ |
| $1/8$ | 9 | 2 | $23/2$ | $2$ |
| $1/10$ | 33 | 8 | $115/2$ | $8$ |
| $2/11$ | 5 | 1 | $27/2$ | $75/77$ |

These are **exact verifications**, not approximate numerical searches.

### 3. Reproduce the 21-term finite rational reduction

```bash
python3 verify_coverage.py
```

**Additional data required:** This script reads `odd_moduli_exact_gaps.csv` and `even_moduli_exact_gaps.csv` from the same directory. **These CSV files were not present in the repository listing when this README was prepared.** Please add the data files before running this command; otherwise it will fail with `FileNotFoundError`.

The expected CSV format is one headerless row per coefficient:

```text
M,p,s,j,gap_n,gap_d
```

representing $G_{p,M}(j/10)=\texttt{gap\\_n}/\texttt{gap\\_d}$. With the data used for the paper, the script constructs admissible coefficient intervals by the Lipschitz bound, checks exact rational interval intersections, and reports the remaining squared scales. The previous research computation gave **5,262** fractions outside the broad inequality and **eight** remaining cases before applying the four explicit colourings above.

After adding the CSV files, independently recheck the qualifying gap values with:

```bash
g++ -O2 -std=c++17 verify_gaps.cpp -o verify_gaps
./verify_gaps
```

`verify_gaps.cpp` includes `gap_core.cpp`; **do not compile both as separate translation units**. It recomputes the gap values using exact arithmetic and checks them against the CSV entries. In the research dataset, this check covered **1,140** qualifying entries.

### 4. Inspect the four remaining 21-term cases

```bash
g++ -O2 -std=c++17 verify_m21_examples.cpp -o verify_m21_examples
./verify_m21_examples
```

This prints exact gap values and a pass/fail comparison for selected colourings at

$$
\alpha^2\in\{\frac14,\frac19,\frac2{15},\frac34\}.
$$

The examples illustrate **why those particular choices fail**. They do **not** prove that no spherical colouring works for the four scales. `verify_m21_examples.cpp` includes `verify_four.cpp`, so compile the former by itself.

### 5. Run the exploratory 21-term search (optional)

```bash
g++ -O2 -std=c++17 m21_target_search.cpp -o m21_target_search
g++ -O2 -std=c++17 m21_verify_screen.cpp -o m21_verify_screen
./m21_target_search > m21_candidates.txt
./m21_verify_screen < m21_candidates.txt
```

The first program searches a restricted family of consecutive-residue spherical colourings and samples the linear coefficient numerically. Its `candidate` lines are **only preliminary leads**. The second program parses these lines and subjects them to an exact collision-based gap test. An absence of successful output does not rule out other moduli, residue sets, or spacings. The search is not exhaustive over all spherical colourings.

## Mathematical background: what is being computed?

For a modulus $p$, integer $M\ge2$, and real $u$, define

$$
G_{p,M}(u) = \max_{\beta \in \mathbb R} \text{gap} \\{ -uk^2 - \beta k \pmod p : 0 \le k < M \\},
$$

where `gap` is the largest circular gap between consecutive points on $\mathbb R/p\mathbb Z$.

The principal colourings are spherical:

$$
x\text{ is red}\quad\Longleftrightarrow\quad
\lfloor R\|x\|^2\rfloor\bmod p\in\{0,1,\ldots,s-1\},
$$

where $R>0$ is selected to avoid red three-term progressions and $U=\alpha^2R$ is the corresponding blue coefficient. A sufficient condition for the blue part is $G_{p,M}(U)\le s$; the red part is checked separately.

The exact computation of $G_{p,M}(u)$ reduces the continuous maximization over $\beta$ to finitely many *collision values*

$$
\beta=\frac{pm}{j-i}-u(i+j),\qquad
0\le i<j<M,\quad 0\le m<j-i.
$$

The code compares rational quantities using integer arithmetic, avoiding floating-point uncertainty in the verification routines. To establish inequalities on whole intervals, the paper combines these exact values with a Lipschitz estimate and finite subdivision.

## File guide

| File | Purpose | Type |
| --- | --- | --- |
| [`verify_uniform22_with_rational17.cpp`](verify_uniform22_with_rational17.cpp) | Main uniform $M=22$ checks and $M=17$ interval checks | Exact verification |
| [`verify_four.cpp`](verify_four.cpp) | Four successful $M=21$ squared-scale examples | Exact verification |
| [`verify_coverage.py`](verify_coverage.py) | Rational-case reduction for $M=21$ | Exact verification; requires CSV data |
| [`verify_gaps.cpp`](verify_gaps.cpp) | Recomputes gap values from the CSV tables | Exact verification; requires CSV data |
| [`gap_core.cpp`](gap_core.cpp) | Shared exact gap routine included by `verify_gaps.cpp` | Supporting source |
| [`verify_m21_examples.cpp`](verify_m21_examples.cpp) | Exact checks for unsuccessful trial colourings | Diagnostic verification |
| [`m21_target_search.cpp`](m21_target_search.cpp) | Samples potential $M=21$ colourings | Exploratory numerical search |
| [`m21_verify_screen.cpp`](m21_verify_screen.cpp) | Exact follow-up test for sampled candidates | Candidate screening |
| [`computing_gpmu.py`](computing_gpmu.py) | Readable Python example computing a circular gap | Educational / exact arithmetic |

For a small standalone example, run:

```bash
python3 computing_gpmu.py
```

The current example computes $G_{29,22}(1/2)$.

## Interpretation and reproducibility

- **Exact verification** means the listed finite arithmetic claims are checked without floating-point approximation. It does not replace the mathematical proof that these checks imply the theorem.
- **Exploratory search** means a restricted numerical investigation. A candidate must pass a separate exact test before it can be used in a proof.
- The scripts using the two CSV tables cannot reproduce the full $M=21$ reduction until those data files are committed alongside the code.
- Some C++ programs use `assert` for verification. Compile **without** `-DNDEBUG`, which would disable those checks.
- If a program appears slow, note that exact gap evaluation enumerates many pairwise coincidences; runtime depends on hardware and compiler optimization.

## Further research

The uniform $M=21$ question remains open within this approach, with the four squared scales $1/4$, $1/9$, $2/15$, and $3/4$ not yet covered by the recorded constructions. Promising directions include exact optimization over the red-admissible spacing intervals and spherical colourings using nonconsecutive residue sets.

## Citation

If you use these programs in research, please cite the accompanying paper (add its full bibliographic information and DOI/arXiv link when available) and link to this repository. The paper contains the mathematical arguments; this repository provides reproducible computations supporting them.

---

**Repository:** https://github.com/Goodbilly/Avoiding-Short-Progressions-Bound
