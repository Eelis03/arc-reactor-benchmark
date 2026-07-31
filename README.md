# Arc Reactor Benchmark

Power balance and confinement benchmarking for ARC-class compact high-field tokamaks with a parametric efficiency study.

[![CI](https://github.com/Eelis03/arc-reactor-benchmark/actions/workflows/ci.yml/badge.svg)](https://github.com/Eelis03/arc-reactor-benchmark/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.12-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

ARC here means the Affordable, Robust, Compact tokamak of Sorbom and colleagues,
published in Fusion Engineering and Design in 2015: a peer-reviewed engineering
design for a compact high-field fusion pilot plant with demountable rare earth
barium copper oxide toroidal field magnets and a molten salt liquid immersion
blanket. It is not the fictional device of a similar nickname, and nothing in
this repository has any relationship to that.

This package solves a zero-dimensional deuterium tritium power balance, evaluates
the operational limits a design point has to respect, converts the result into
net electrical output, and compares all of it against three published design
points: ARC, ITER, and SPARC. It is written for someone who wants to know how far
a compact high-field concept can be pushed before something binds, and which
number in the projection is carrying the most uncertainty.

## Problem

A fusion power plant projection is a chain of multiplications, and the answer at
the end of it is far more sensitive to the inputs than the inputs are to each
other. The energy confinement time comes from an empirical scaling fitted to
existing machines and extrapolated well outside them. That confinement time sets
the loss power, the loss power sets the required auxiliary heating, and the
auxiliary heating is drawn from the grid at less than sixty percent efficiency.

The amplification is measurable and it is large. At a solved steady state the
loss power goes as the stored energy raised to `1 / (1 - alpha_P)`, with
`alpha_P` the power degradation exponent of the confinement scaling. For
IPB98(y,2) that exponent is 0.69, so the amplification factor is 3.23: a ten
percent error in the assumed confinement enhancement moves the loss power by
roughly a third, and moves the fusion gain by more. A projection that reports a
gain without reporting what it assumed about confinement, and how sensitive the
answer is to that assumption, is reporting the assumption rather than a result.

The second problem is that a plasma result is not a plant result. A device can
reach a respectable fusion gain and still consume more electricity than it
produces, because the auxiliary heating that holds the plasma is bought at the
wall plug and sold back at the thermal conversion efficiency. This package
computes both, and reports the case where they disagree.

## Approach

The plasma is described by one volume-averaged density, one volume-averaged
temperature, a composition, and a shape. Fusion power uses the Bosch and Hale
(1992) Maxwellian reactivity fit for T(d,n)4He. Losses are bremsstrahlung with
the relativistic and electron-electron correction of Rider (1995), synchrotron
radiation as the Larmor emission multiplied by a Trubnikov escape factor,
impurity line radiation from caller-supplied coronal loss parameters, and
transport through the energy confinement time.

The confinement time comes from a scaling law reached through a Protocol, so the
solver never names a particular fit. Three are provided: IPB98(y,2) from the ITER
Physics Basis, the ITER89-P L-mode scaling of Yushmanov and colleagues, and the
constant-beta H-mode fit of Petty (2008). Every conclusion in the results section
is reported against all three.

Restricting the Protocol to scalings that are power laws in the loss power buys
something specific. The steady-state balance then reduces to
`P_loss = (W / tau_1)**(1 / (1 - alpha_P))` and is solved in closed form, with no
iteration, no tolerance, and no convergence test. The only iterative solves in
the package are a bracketed root find for an ignition temperature and a bounded
minimisation for the Lawson optimum, and both verify convergence before
returning.

Four operational limits are evaluated as constraints rather than as commentary:
the Greenwald density limit, the Troyon beta limit, the safety factor at the 95
percent flux surface, and a bootstrap current fraction estimate. A fifth, the L to
H transition power threshold of Martin and colleagues (2008), is applied whenever
an H-mode scaling is used, because an H-mode confinement time is not available to
a plasma that cannot reach H-mode. A design point that violates any of them is
reported as violating it. The alternatives that were considered and rejected are
recorded in [docs/design-notes.md](docs/design-notes.md).

## Installation

Requires Python 3.12 or later.

```bash
git clone https://github.com/Eelis03/arc-reactor-benchmark.git
cd arc-reactor-benchmark
uv sync
```

Using pip instead of uv:

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

## Usage

```python
from arc_benchmark import IPB98Y2, evaluate_constraints, machine, solve_operating_point

case = machine("ARC")
point = solve_operating_point(case.state, IPB98Y2)
report = evaluate_constraints(point)

print(f"fusion power   {point.terms.fusion_power_mw:.1f} MW")
print(f"auxiliary      {point.auxiliary_power_mw:.1f} MW")
print(f"fusion gain Q  {point.fusion_gain:.2f}")
print(f"binding limit  {report.binding.name} at {report.binding.utilisation:.3f}")
```

Swapping `IPB98Y2` for `PETTY08` or `ITER89P` is the only change needed to see the
same design point under a different confinement scaling, and swapping `"ARC"` for
`"ITER"` or `"SPARC"` is the only change needed to move to another machine.

Runnable examples live in `examples/`:

```bash
uv run python examples/power_balance.py
uv run python examples/benchmark_points.py
uv run python examples/field_scaling.py
uv run python examples/efficiency_study.py
uv run python examples/scaling_sensitivity.py
```

Each accepts `--quick` for a short run without figures, and the sweep scripts
accept `--points` to change the resolution.

## Results

Every number below is the output of the command named above it. The default
configuration throughout is flat profiles, the IPB98(y,2) scaling, the separatrix
loss power convention, and a first wall reflectivity of 0.9. Flat profiles are
the honest zero-dimensional answer, and the size of what they leave out is
reported separately rather than absorbed into a fitted factor.

### Benchmark against published design points

From `uv run python examples/benchmark_points.py`. The published values are design
targets from the cited papers, not measurements: none of the three machines has
operated at these points.

| Machine | Quantity | Computed | Published | Error |
| --- | --- | --- | --- | --- |
| ARC | fusion power, MW | 355.1 | 525.0 | -32.4 percent |
| ARC | fusion gain Q | 3.16 | 13.6 | -76.7 percent |
| ARC | auxiliary power, MW | 112.2 | 38.6 | +190.7 percent |
| ARC | normalised beta | 2.250 | 2.590 | -13.1 percent |
| ARC | bootstrap fraction | 0.647 | 0.630 | +2.7 percent |
| ARC | Greenwald fraction | 0.669 | 0.670 | -0.2 percent |
| ARC | thermal power, MW | 552.6 | 708.0 | -22.0 percent |
| ARC | net electric power, MWe | -6.9 | 190.0 | -103.6 percent |
| ITER | fusion power, MW | 335.6 | 500.0 | -32.9 percent |
| ITER | fusion gain Q | 3.66 | 10.0 | -63.4 percent |
| ITER | confinement time, s | 3.090 | 3.700 | -16.5 percent |
| ITER | auxiliary power, MW | 91.7 | 50.0 | +83.4 percent |
| ITER | normalised beta | 1.693 | 1.800 | -5.9 percent |
| ITER | bootstrap fraction | 0.224 | 0.150 | +49.4 percent |
| ITER | Greenwald fraction | 0.838 | 0.850 | -1.4 percent |
| SPARC | fusion power, MW | 66.4 | 140.0 | -52.6 percent |
| SPARC | fusion gain Q | 1.84 | 11.0 | -83.3 percent |
| SPARC | confinement time, s | 0.788 | 0.770 | +2.3 percent |
| SPARC | auxiliary power, MW | 36.2 | 11.0 | +229.0 percent |
| SPARC | Greenwald fraction | 0.364 | 0.370 | -1.7 percent |

The agreements and the disagreements separate cleanly, and the pattern is more
informative than any single row.

Quantities that are closed forms of the machine parameters are reproduced to
within two percent: the Greenwald fraction is out by 0.2, 1.4, and 1.7 percent
for the three machines. Quantities that depend on the pressure are out by five to
thirteen percent, which is the size of the profile and fast alpha contributions
this model does not carry. The bootstrap fraction agrees for ARC to 2.7 percent
and is 49 percent high for ITER, which is what a single coefficient of 0.7 in
`f_BS = c sqrt(eps) beta_p` does when applied across a factor of seven in
poloidal beta.

The fusion power is 32 to 53 percent low everywhere, and that is the flat profile
assumption doing exactly what it should: `<n**2 sigma_v(T)>` exceeds
`<n>**2 sigma_v(<T>)` for any peaked profile, and this model evaluates the second.

The fusion gain is 63 to 83 percent low, a larger error than the fusion power
alone accounts for. The reason is the amplification described in the problem
statement: the auxiliary power is the difference between the transport loss and
the alpha heating, so an error in the alpha heating and an error in the loss power
both land on it, and the loss power carries the exponent 3.23.

SPARC as published is reported as violating a constraint. Its solved loss power
of 31.8 MW sits below the 40.7 MW L to H access threshold, a shortfall of a
factor of 1.28, so the H-mode confinement time it was solved with is not
available to it under this model's own rules. That is recorded as a violation
rather than a footnote, and it is pinned as such in the regression tests.

### The implied confinement enhancement, which separates the two failures

The fusion gain conflates a fusion power disagreement with a confinement
disagreement. Taking the published fusion power and the published auxiliary power
as given, and asking what H factor this model would need to reproduce them,
separates the two.

| Machine | H factor assumed by the source | H factor implied by the published point | Ratio |
| --- | --- | --- | --- |
| ARC | 1.800 | 2.035 | 1.130 |
| ITER | 1.000 | 1.027 | 1.027 |
| SPARC | 1.000 | 1.132 | 1.132 |

The confinement model agrees with all three sources to within 13 percent, and
with ITER to 2.7 percent. ITER is the point IPB98(y,2) was constructed to project,
so that agreement is the check that the scaling is transcribed correctly rather
than a claim about the model. Evaluated directly at the ITER reference parameters
with a loss power of 87 MW, IPB98(y,2) here returns 3.59 s against the published
3.7 s, a difference of 2.9 percent; the published loss power is quoted between 80
and 90 MW across sources, which alone moves the result from 3.81 s to 3.51 s.

The conclusion is that the gain discrepancies are a fusion power problem, not a
confinement problem. That is worth knowing, because the fusion power problem has a
known cause and the confinement one would not have.

### What profiles are worth, and why one shape does not serve three machines

The same benchmark with parabolic profiles at a density exponent of 0.4 and a
temperature exponent of 1.0, which enhances fusion power by 1.65 and stored energy
by 1.17 at 14 keV:

| Machine | Flat, MW | Peaked, MW | Published, MW |
| --- | --- | --- | --- |
| ARC | 355.1 | 587.4 | 525.0 |
| ITER | 335.6 | 771.8 | 500.0 |
| SPARC | 66.4 | 176.1 | 140.0 |

One profile shape moves ARC from 32 percent low to 12 percent high, ITER from 33
percent low to 54 percent high, and SPARC from 53 percent low to 26 percent high.
ARC's normalised beta improves from 13 percent low to 1.4 percent high in the same
step. The shape that suits one machine overshoots another, which is the direct
statement that profile shape is a per-machine input and not a universal
correction. Choosing a shape per machine to match the published fusion power would
turn the benchmark into a fit, so the flat case is what the primary table reports.

### Field scaling, quantified

From `uv run python examples/field_scaling.py`, 41 points from 4 T to 14 T at
fixed size and shape, with the plasma current tracked so that the cylindrical
safety factor stays at 5.004.

| Held fixed | Fusion power exponent in field | Analytic expectation | Gain exponent | Confinement time exponent |
| --- | --- | --- | --- | --- |
| Electron density | +0.0000 | 0 | +2.5024 | +3.4839 |
| Greenwald fraction | +2.0000 | 2 | +3.2664 | +2.5806 |
| Toroidal beta | +4.0000 | 4 | +4.4189 | +1.6774 |

All three fits return a coefficient of determination of 1.000000, because the
underlying relation is an exact power law in each case. The exponents are the
whole high-field argument stated as a number. At fixed density the field buys
nothing directly in fusion power and everything through confinement. At a fixed
Greenwald fraction the density tracks the current, which tracks the field, so
fusion power goes as the square of the field. At fixed beta the density goes as
the square of the field and fusion power as the fourth power, which is the branch
the compact high-field case is normally argued on.

The gain exponents exceed the fusion power exponents in every case, because
raising the field also lengthens the confinement time and so reduces the heating
required.

Which limit binds first depends on the direction:

| Held fixed | Feasible points | Binding constraint tally | Where it fails |
| --- | --- | --- | --- |
| Electron density | 12 of 41, from 8.25 T to 11.00 T | Troyon at 25, L to H threshold at 16 | Troyon, Greenwald, and bootstrap below 8.25 T; L to H threshold above 11.00 T |
| Greenwald fraction | 18 of 41, from 7.50 T to 11.75 T | Troyon at 27, L to H threshold at 14 | Troyon and bootstrap below 7.50 T; L to H threshold above 11.75 T |
| Toroidal beta | 40 of 41, from 4.00 T to 13.75 T | Troyon at 29, Greenwald at 12 | Greenwald and the L to H threshold at 14.00 T |

The fixed-beta branch is not only the steepest, it is also the one that stays
feasible longest, because holding beta constant is precisely holding the Troyon
constraint at a fixed utilisation while the field rises. Its failure at 14 T is
the Greenwald limit: the density needed to hold beta at that field exceeds what
the current density supports.

Sweeping the major radius instead, from 2 m to 8 m at fixed aspect ratio, field,
and safety factor, gives a fusion power exponent of +1.0000 against the analytic
expectation of 1. The volume goes as the cube of the radius while the Greenwald
density falls as its inverse, and fusion power carries the square of the density.
Only 12 of 41 points are feasible, from 2.75 m to 4.40 m: Troyon and the bootstrap
fraction fail below that band and the L to H access threshold fails above it. That
threshold is the binding constraint at 28 of the 41 points, because a large device
at a fixed Greenwald fraction has a long confinement time, therefore a low loss
power, and therefore cannot reach H-mode.

### Net electrical efficiency

From `uv run python examples/efficiency_study.py`, at the flat-profile ARC point
with a thermal efficiency of 0.40, a blanket energy multiplication of 1.30, a
heating wall plug efficiency of 0.55, a 10 MWe cryoplant, a balance of plant at 4
percent of gross, and 5 MWe of tritium and site load.

| Term | Value |
| --- | --- |
| Fusion power | 355.10 MW |
| Neutron power | 284.08 MW |
| Thermal power | 552.55 MW |
| Gross electric | 221.02 MWe |
| Heating and current drive draw | 204.04 MWe |
| Cryoplant | 10.00 MWe |
| Balance of plant | 8.84 MWe |
| Tritium plant and site | 5.00 MWe |
| Recirculating total | 227.88 MWe |
| Net electric | -6.86 MWe |
| Recirculating fraction | 1.031 |
| Engineering gain | 0.970 |

The flat-profile ARC point does not produce net electricity in this model, and it
is not close to producing it: the heating and current drive draw alone, 204 MWe,
is 92 percent of the gross output. That number is the whole efficiency problem in
one line, and it follows directly from the fusion gain of 3.16.

Where the sensitivity lies, one parameter at a time from the same baseline:

| Parameter | Value | Net electric, MWe | Net efficiency |
| --- | --- | --- | --- |
| Thermal efficiency | 0.33 | -43.99 | -0.1239 |
| Thermal efficiency | 0.40 | -6.86 | -0.0193 |
| Thermal efficiency | 0.45 | 19.66 | 0.0554 |
| Thermal efficiency | 0.50 | 46.18 | 0.1301 |
| Heating wall plug efficiency | 0.35 | -123.46 | -0.3477 |
| Heating wall plug efficiency | 0.45 | -52.21 | -0.1470 |
| Heating wall plug efficiency | 0.55 | -6.86 | -0.0193 |
| Heating wall plug efficiency | 0.70 | 36.86 | 0.1038 |
| Blanket multiplication | 1.10 | -28.68 | -0.0808 |
| Blanket multiplication | 1.40 | 4.05 | 0.0114 |
| Cryoplant | 2 MWe | 1.14 | 0.0032 |
| Cryoplant | 50 MWe | -46.86 | -0.1320 |

Moving the heating wall plug efficiency from 0.35 to 0.70 is worth 160 MWe, more
than any other parameter in the table. That is a consequence of the low gain, and
it is the reason a current drive system is the component a recirculating-power
study should look at first at this operating point. The cryoplant, often raised as
the objection to a superconducting machine, is worth 48 MWe across the whole range
from 2 to 50 MWe, a factor of three less than the heating system.

The first wall reflectivity acts on the plasma rather than on the plant, and it
acts strongly, because the synchrotron loss carries `sqrt(1 - R_w)`:

| Reflectivity | Synchrotron, MW | Auxiliary, MW | Q | Net electric, MWe |
| --- | --- | --- | --- | --- |
| 0.60 | 111.59 | 168.02 | 2.11 | -86.88 |
| 0.80 | 78.91 | 135.33 | 2.62 | -40.01 |
| 0.90 | 55.80 | 112.22 | 3.16 | -6.86 |
| 0.95 | 39.45 | 95.88 | 3.70 | 16.57 |
| 0.98 | 24.95 | 81.38 | 4.36 | 37.37 |

At the baseline reflectivity the synchrotron loss of 55.80 MW is nine times the
bremsstrahlung loss of 6.06 MW and is by far the largest radiated term. It is also
the least reliable number in the model, for the reasons set out in the design
notes, so the sensitivity is reported rather than a single value. Whether this
design point produces net electricity at all depends on it.

Sweeping the density from 0.6e20 to 2.4e20 per cubic metre with everything else
held gives a best feasible net efficiency of 0.0629 at 1.59e20 per cubic metre,
with 21 of 41 points feasible, spanning 0.69e20 to 1.59e20 per cubic metre. The
unconstrained optimum is 0.1743 at 2.40e20 per cubic metre, which violates the
Troyon, Greenwald, and bootstrap limits and is reported as violating them. Below
0.69e20 per cubic metre the loss power falls under the L to H access threshold, so
the feasible band is bounded on both sides.

### Sensitivity to the confinement scaling

From `uv run python examples/scaling_sensitivity.py`.

| Scaling | Mode | Power degradation | Amplification into loss power |
| --- | --- | --- | --- |
| IPB98(y,2) | H | 0.690 | 3.226 |
| ITER89-P | L | 0.500 | 2.000 |
| Petty08 | H | 0.470 | 1.887 |

The same ARC design point under each:

| Scaling | Confinement time, s | Auxiliary power, MW | Q | Implied H |
| --- | --- | --- | --- | --- |
| IPB98(y,2) | 1.075 | 112.22 | 3.16 | 2.035 |
| ITER89-P | 0.418 | 303.10 | 1.17 | 3.518 |
| Petty08 | 2.571 | 41.61 | 8.53 | 1.399 |

The gain moves by a factor of 7.3 across three published scalings applied to one
unchanged design point. This is the single largest source of uncertainty in the
projection, larger than the profile shape and far larger than any engineering
parameter in the efficiency table.

Holding the scaling fixed and varying only the assumed H factor:

| H factor | Confinement time, s | Loss power, MW | Auxiliary power, MW | Q |
| --- | --- | --- | --- | --- |
| 1.20 | 0.291 | 448.98 | 439.81 | 0.81 |
| 1.50 | 0.597 | 218.58 | 209.41 | 1.70 |
| 1.80 | 1.075 | 121.39 | 112.22 | 3.16 |
| 2.10 | 1.768 | 73.83 | 64.66 | 5.49 |
| 2.50 | 3.103 | 42.07 | 32.90 | 10.79 |
| 3.00 | 5.587 | 23.36 | 14.20 | 25.01 |

Raising the H factor from 1.8 to 2.5, a 39 percent change, raises the gain from
3.16 to 10.79, a factor of 3.4. The measured ratio of loss powers, 121.39 to
42.07, is 2.885, against the analytic prediction of `(2.5 / 1.8)**3.226`, which is
2.885. The amplification exponent is not an estimate.

The choice of loss power convention is worth almost as much as a change of
scaling. Evaluating IPB98(y,2) at the power crossing the separatrix, which
subtracts core radiation, gives an auxiliary power of 112.22 MW and a gain of
3.16. Evaluating it at the total heating power, which is the definition the
IPB98 regression itself was fitted with, gives 50.37 MW and a gain of 7.05. Both
are defensible and the difference is a factor of 2.2 in the reported gain, so the
convention is a parameter here rather than a hidden decision.

### Lawson condition

From `uv run python examples/power_balance.py`. Derived from the same balance the
rest of the package solves, rather than quoted.

For a pure deuterium tritium plasma with classical bremsstrahlung and no
dilution, the minimum ignition triple product is 2.95e21 per cubic metre keV
second at 14.25 keV, against the roughly 3e21 quoted in the standard references.
Adding 5 percent helium ash and the relativistic bremsstrahlung correction raises
it to 3.65e21 at 14.50 keV, a 24 percent penalty for a 10 percent fuel dilution,
because fusion power carries the square of the fuel fraction while the stored
energy carries only the first power.

The flat-profile ARC point achieves a triple product of 1.957e21 per cubic metre
keV second, which is 0.535 of the ignition requirement at its own temperature and
composition. It is reported as not ignited, and the ignition temperature search
returns no root, because with the IPB98(y,2) power degradation the transport loss
rises as the temperature to the power 3.23 while the reactivity flattens above 30
keV. At this density, ignition is not something raising the temperature can reach.

Synchrotron radiation is deliberately excluded from the Lawson condition. It
carries a half power of density rather than the square, so including it would
make the condition a function of density and the curve would stop being a curve.
Bremsstrahlung and coronal line radiation both scale as the square of density and
are included.

## Architecture

| Module | Responsibility |
| --- | --- |
| `src/arc_benchmark/_types.py` | The float64 array alias every public signature uses |
| `src/arc_benchmark/model/constants.py` | Physical constants and the deuterium tritium reaction energetics |
| `src/arc_benchmark/model/geometry.py` | Plasma shape: volume, surface area, poloidal perimeter, aspect ratio |
| `src/arc_benchmark/model/reactivity.py` | Bosch and Hale reactivity fit, cross-section fit, and the Maxwellian integral of the second |
| `src/arc_benchmark/model/profiles.py` | Parabolic profile correction factors, closed form except for fusion power |
| `src/arc_benchmark/model/radiation.py` | Bremsstrahlung, cyclotron emission with the Trubnikov escape factor, impurity line radiation |
| `src/arc_benchmark/model/confinement.py` | IPB98(y,2), ITER89-P, and Petty08 in one canonical exponent basis |
| `src/arc_benchmark/model/limits.py` | Greenwald, Troyon, safety factor, bootstrap fraction, L to H threshold |
| `src/arc_benchmark/model/plant.py` | Thermal conversion, recirculating power, net efficiency, engineering gain |
| `src/arc_benchmark/algorithm/protocols.py` | The `ConfinementScaling` Protocol the solver works against |
| `src/arc_benchmark/algorithm/balance.py` | Plasma state, composition, and every state-dependent power term |
| `src/arc_benchmark/algorithm/operating.py` | Closed-form steady-state solve and the bracketed ignition temperature search |
| `src/arc_benchmark/algorithm/constraints.py` | Limits turned into verdicts, with utilisation and the binding constraint |
| `src/arc_benchmark/algorithm/lawson.py` | Ignition condition, gain-dependent requirement, and the triple product optimum |
| `src/arc_benchmark/pipeline/machines.py` | ARC, ITER, and SPARC with their published values and their sources |
| `src/arc_benchmark/pipeline/sweep.py` | Field, radius, and density sweeps with an explicit invariant |
| `src/arc_benchmark/pipeline/benchmark.py` | Reproduction of published points and the implied H factor |
| `src/arc_benchmark/pipeline/trace.py` | The structured records the pipeline produces |
| `src/arc_benchmark/analysis/metrics.py` | Log-space power law fits and binding constraint tallies |
| `src/arc_benchmark/analysis/report.py` | Text tables built from traces |
| `src/arc_benchmark/analysis/tables.py` | The same traces as data frames |
| `src/arc_benchmark/analysis/figures.py` | Sweep, Lawson, and power balance figures |
| `examples/` | Thin wiring scripts with no logic of their own |

The dependency direction is one way: `model` knows nothing of anything else,
`algorithm` uses `model`, `pipeline` uses both, `analysis` reads what `pipeline`
produces, and `examples` only wires those together. The confinement scalings live
in `model` and satisfy the Protocol declared in `algorithm` structurally, without
importing it, so a caller can supply a scaling this package has never heard of.

## Testing

```bash
uv run pytest
uv run ruff check .
uv run mypy
```

The suite has three tiers: property and invariant tests covering the mathematics,
regression tests pinning recorded behaviour, and integration tests running each
example script under a reduced iteration count.

Tier one covers the mathematics, with four checks that reach outside this package
entirely. The Bosch and Hale reactivity is compared against the values tabulated
in the source paper at five temperatures spanning four orders of magnitude. The
same reactivity is compared against a numerical Maxwellian average of the separate
cross-section fit from the same paper, which has different coefficients and was
derived by a different procedure; the two agree to within 0.75 percent across 1 to
100 keV. IPB98(y,2) is compared against the published ITER confinement time. The
Greenwald fraction, the safety factor at the 95 percent flux surface, and the L to
H access threshold are compared against published values for ITER, ARC, and SPARC.
The coefficient of the cyclotron emission is assembled from the elementary charge,
the electron mass, the electric constant, and the speed of light, and checked
against the literal usually quoted for it.

The rest of tier one asserts properties by construction: the power balance closes
to 1e-9 MW under both loss power conventions and all three scalings; the solved
confinement time is the one the scaling returns at the solved loss power; the gain
rises monotonically with the confinement time; the ignition boundary, located by
root finding, has the alpha power exactly equal to radiation plus transport;
doubling the field at fixed beta multiplies fusion power by exactly sixteen and at
a fixed Greenwald fraction by exactly four, both against the analytic expectation;
every confinement exponent is recovered by perturbing its own input; and a design
point above the Greenwald limit or below the safety factor limit is flagged rather
than returned.

Tier two pins a recorded reference run. Every operating point in this package is a
closed-form inversion with no iteration, so those values are pinned at 1e-10
relative, a tolerance derived from the operation count and the machine epsilon.
The one iterative result that is pinned, the Lawson optimum, is pinned to the
tolerance its minimiser was asked to converge to and not to the digits it happens
to produce. Sweep-derived quantities are quantised by the sweep resolution, so a
crossing field is pinned to one sweep step and counts are pinned as integers. An
unconverged solve is never pinned: the ignition temperature search returns `None`
when its bracket contains no sign change, and that is what the test asserts. The
reasoning is written out in the module docstring of `tests/test_regression.py`.

Tier three loads each script in `examples/` and calls its `main` with `--quick`,
asserting a zero return code, non-empty output, and the specific facts each script
exists to produce. One script is additionally run with figures enabled so that the
plotting code is exercised.

The full suite of 206 tests runs in under 2 seconds.

## References

Physics:

- Sorbom, B. N., J. Ball, T. R. Palmer, F. J. Mangiarotti, J. M. Sierchio,
  P. Bonoli, C. Kasten, D. A. Sutherland, H. S. Barnard, C. B. Haakonsen,
  J. Goh, C. Sung, and D. G. Whyte. "ARC: A compact, high-field, fusion nuclear
  science facility and demonstration power plant with demountable magnets."
  *Fusion Engineering and Design* 100 (2015): 378 to 405.
  DOI [10.1016/j.fusengdes.2015.07.008](https://doi.org/10.1016/j.fusengdes.2015.07.008).
  The design point this repository is named after and benchmarked against.
- Bosch, H.-S., and G. M. Hale. "Improved formulas for fusion cross-sections and
  thermal reactivities." *Nuclear Fusion* 32, no. 4 (1992): 611 to 631.
  DOI [10.1088/0029-5515/32/4/I07](https://doi.org/10.1088/0029-5515/32/4/I07).
  Both the T(d,n)4He cross-section parameterisation and the Maxwellian reactivity
  parameterisation, and the tabulated values the test suite checks against.
- ITER Physics Expert Groups on Confinement and Transport and Confinement
  Modelling and Database. "Chapter 2: Plasma confinement and transport."
  *Nuclear Fusion* 39, no. 12 (1999): 2175 to 2249.
  DOI [10.1088/0029-5515/39/12/302](https://doi.org/10.1088/0029-5515/39/12/302).
  The IPB98(y,2) ELMy H-mode thermal energy confinement scaling and the bootstrap
  fraction discussion.
- ITER Physics Basis Editors. "Chapter 1: Overview and summary."
  *Nuclear Fusion* 39, no. 12 (1999): 2137 to 2174.
  DOI [10.1088/0029-5515/39/12/301](https://doi.org/10.1088/0029-5515/39/12/301).
  The inductive Q equal to ten reference scenario used as the ITER benchmark
  point.
- Greenwald, M., J. L. Terry, S. M. Wolfe, S. Ejima, M. G. Bell, S. M. Kaye, and
  G. H. Neilson. "A new look at density limits in tokamaks." *Nuclear Fusion* 28,
  no. 12 (1988): 2199 to 2207.
  DOI [10.1088/0029-5515/28/12/009](https://doi.org/10.1088/0029-5515/28/12/009).
  The density limit.
- Troyon, F., R. Gruber, H. Saurenmann, S. Semenzato, and S. Succi. "MHD limits to
  plasma confinement." *Plasma Physics and Controlled Fusion* 26, no. 1A (1984):
  209 to 215.
  DOI [10.1088/0741-3335/26/1A/319](https://doi.org/10.1088/0741-3335/26/1A/319).
  The normalised beta limit.
- Lawson, J. D. "Some criteria for a power producing thermonuclear reactor."
  *Proceedings of the Physical Society. Section B* 70, no. 1 (1957): 6 to 10.
  DOI [10.1088/0370-1301/70/1/303](https://doi.org/10.1088/0370-1301/70/1/303).
  The condition the triple product analysis derives from.
- Yushmanov, P. N., T. Takizuka, K. S. Riedel, O. J. W. F. Kardaun,
  J. G. Cordey, S. M. Kaye, and D. E. Post. "Scalings for tokamak energy
  confinement." *Nuclear Fusion* 30, no. 10 (1990): 1999 to 2006.
  DOI [10.1088/0029-5515/30/10/001](https://doi.org/10.1088/0029-5515/30/10/001).
  The ITER89-P L-mode scaling, used here as the pessimistic bound.
- Petty, C. C. "Sizing up plasmas using dimensionless parameters."
  *Physics of Plasmas* 15, no. 8 (2008): 080501.
  DOI [10.1063/1.2961043](https://doi.org/10.1063/1.2961043).
  The constant-beta constrained H-mode confinement fit.
- Martin, Y. R., T. Takizuka, and the ITPA CDBM H-mode Threshold Database Working
  Group. "Power requirement for accessing the H-mode in ITER."
  *Journal of Physics: Conference Series* 123 (2008): 012033.
  DOI [10.1088/1742-6596/123/1/012033](https://doi.org/10.1088/1742-6596/123/1/012033).
  The L to H transition power threshold applied as a constraint on every H-mode
  solve.
- Trubnikov, B. A. "Universal coefficients for synchrotron emission from plasma
  configurations." In *Reviews of Plasma Physics*, volume 7, edited by
  M. A. Leontovich, 345 to 379. New York: Consultants Bureau, 1979.
  DOI [10.1007/978-1-4757-0836-4](https://doi.org/10.1007/978-1-4757-0836-4).
  The asymptotic cyclotron escape factor used for the synchrotron term.
- Albajar, F., J. Johner, and G. Granata. "Improved calculation of synchrotron
  radiation losses in realistic tokamak plasmas." *Nuclear Fusion* 41, no. 6
  (2001): 665 to 678.
  DOI [10.1088/0029-5515/41/6/301](https://doi.org/10.1088/0029-5515/41/6/301).
  Cited as the accurate treatment that this package does not implement, and as
  the reason the Trubnikov form is expected to overestimate.
- Rider, T. H. "A general critique of inertial-electrostatic confinement fusion
  systems." *Physics of Plasmas* 2, no. 6 (1995): 1853 to 1872.
  DOI [10.1063/1.871273](https://doi.org/10.1063/1.871273).
  The relativistic and electron-electron corrections to bremsstrahlung.
- Wilson, H. R. "Bootstrap current scaling in tokamaks." *Nuclear Fusion* 32,
  no. 2 (1992): 257 to 263.
  DOI [10.1088/0029-5515/32/2/I05](https://doi.org/10.1088/0029-5515/32/2/I05).
  The origin of the bootstrap fraction estimate used here.
- Sauter, O., C. Angioni, and Y. R. Lin-Liu. "Neoclassical conductivity and
  bootstrap current formulas for general axisymmetric equilibria and arbitrary
  collisionality regime." *Physics of Plasmas* 6, no. 7 (1999): 2834 to 2839.
  DOI [10.1063/1.873240](https://doi.org/10.1063/1.873240).
  Cited as the neoclassical bootstrap calculation that this package does not
  implement, and the reason it does not.
- Uckan, N. A., and the ITER Physics Group. *ITER Physics Design Guidelines: 1989.*
  ITER Documentation Series No. 10. Vienna: IAEA, 1990.
  [Stable URL](https://inis.iaea.org/records/8jasy-vpp88).
  The shaping-corrected formula for the safety factor at the 95 percent flux
  surface.
- Wesson, J. *Tokamaks.* 4th edition. Oxford: Oxford University Press, 2011.
  ISBN 978-0-19-959223-4. The classical bremsstrahlung coefficient and the
  cylindrical safety factor.
- Creely, A. J., M. J. Greenwald, S. B. Ballinger, D. Brunner, J. Canik,
  J. Doody, T. Fuelop, D. T. Garnier, R. Granetz, T. K. Gray, C. Holland,
  N. T. Howard, J. W. Hughes, J. H. Irby, V. A. Izzo, G. J. Kramer, A. Q. Kuang,
  B. LaBombard, Y. Lin, B. Lipschultz, N. C. Logan, J. D. Lore, E. S. Marmar,
  K. Montes, R. T. Mumgaard, C. Paz-Soldan, C. Rea, M. L. Reinke,
  P. Rodriguez-Fernandez, K. Saerkimaeki, F. Sciortino, S. D. Scott, A. Snicker,
  P. B. Snyder, B. N. Sorbom, R. Sweeney, R. A. Tinguely, E. A. Tolman,
  M. Umansky, O. Vallhagen, J. Varje, D. G. Whyte, J. C. Wright, S. J. Wukitch,
  and J. Zhu. "Overview of the SPARC tokamak." *Journal of Plasma Physics* 86,
  no. 5 (2020): 865860502.
  DOI [10.1017/S0022377820001257](https://doi.org/10.1017/S0022377820001257).
  The third benchmark point.
- Puetterich, T., R. Neu, R. Dux, A. D. Whiteford, M. G. O'Mullane,
  H. P. Summers, and the ASDEX Upgrade Team. "Calculation and experimental test
  of the cooling factor of tungsten." *Nuclear Fusion* 50, no. 2 (2010): 025012.
  DOI [10.1088/0029-5515/50/2/025012](https://doi.org/10.1088/0029-5515/50/2/025012).
  The source of the order-of-magnitude coronal loss parameters supplied as
  defaults, which are not a substitute for an atomic database.

Dependencies:

- [NumPy](https://numpy.org/) (>= 2.0), BSD 3-Clause. Array arithmetic for the
  vectorised reactivity and radiation evaluations, and the log-space least
  squares fit behind every measured exponent.
- [SciPy](https://scipy.org/) (>= 1.14), BSD 3-Clause. Three routines only: the
  composite Simpson rule for the Maxwellian cross-section integral and the
  profile factors, Brent root finding for the ignition temperature, and bounded
  scalar minimisation for the Lawson optimum.
- [pandas](https://pandas.pydata.org/) (>= 2.2), BSD 3-Clause. Renders sweep,
  constraint, and benchmark traces as data frames so that a study can sort,
  filter, and join them without walking the trace objects.
- [Matplotlib](https://matplotlib.org/) (>= 3.9), Matplotlib license, a BSD-style
  permissive license. Sweep, Lawson, and power balance figures, used under the
  non-interactive Agg backend.
- [pytest](https://pytest.org/) (>= 8.3), MIT. Test runner, development only.
- [Ruff](https://docs.astral.sh/ruff/) (>= 0.8), MIT. Linter and import sorter,
  development only.
- [mypy](https://mypy-lang.org/) (>= 1.13), MIT. Static type checker in strict
  mode, development only.

## License

Released under the MIT license. See [LICENSE](LICENSE).
