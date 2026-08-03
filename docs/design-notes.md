# Design notes for Arc Reactor Benchmark

## Method selection

### A zero-dimensional balance, and what that decision costs

The model solves one volume-averaged power balance,

```
P_alpha + P_aux = P_bremsstrahlung + P_synchrotron + P_line + W / tau_E
```

with the plasma described by a single density, a single temperature shared by
electrons and ions, a composition, and a shape. There is no radial coordinate.

That is a scope decision, not a simplification of something larger that was
attempted first. A zero-dimensional balance is the right instrument for the
question this repository asks, which is how the answer moves when the field, the
size, the density, or the confinement assumption moves. It is the wrong
instrument for the question of what the profiles are, and it does not attempt
that. The section "Known limitations" below states exactly what this excludes.

### The steady state is solved in closed form, deliberately

The confinement scaling has the form `tau_E = tau_1 P_loss**(-alpha_P)`. In
steady state the transport loss equals the stored energy divided by the
confinement time, and the confinement time is evaluated at that same loss power,
which gives

```
P_loss = W / (tau_1 P_loss**(-alpha_P))   ==>   P_loss = (W / tau_1)**(1 / (1 - alpha_P))
```

There is no fixed point iteration, no convergence tolerance, and no stopping
test. This was chosen over a general nonlinear solve for two reasons.

The first is that it makes the amplification visible as an exponent rather than
hidden inside an iteration. The derivative `d ln P_loss / d ln tau_1` is exactly
`-1 / (1 - alpha_P)`, which is 3.226 for IPB98(y,2). The results section reports
that as a measured ratio, 121.39 MW against 42.07 MW when the H factor moves from
1.8 to 2.5, against an analytic 2.885. A general solver would produce the same
numbers and would not make the structure legible.

The second reason is reproducibility, and it is the reason the interface is
restricted rather than merely the implementation. A converged iterative solve is
reproducible; a solve that has not converged is not, because after enough
iterations a difference in the order a linear algebra kernel reduces a sum grows
into a visibly different answer. Forcing the confinement scaling through a
Protocol that exposes its power exponent means every operating point in this
package is a handful of arithmetic operations deep, and the regression tests can
pin them at 1e-10 relative without the pin becoming a statement about the machine
that recorded it.

The price is that the Protocol excludes confinement models that are not power
laws in the loss power. A neural network surrogate, or a scaling with a
saturation term, would need a root find added. That is a real restriction and it
is stated in the Protocol docstring rather than discovered later.

### Three confinement scalings, all reported

IPB98(y,2) is the required scaling and is implemented exactly as the ITER Physics
Basis publishes it, in the published units, with no translation. ITER89-P and
Petty08 are provided alongside because the choice between them changes the answer
by more than anything else in the model: the ARC gain moves from 1.17 to 3.16 to
8.53 across the three, a factor of 7.3, on an unchanged design point.

They are written in one canonical exponent basis so that the same solver and the
same tests apply to all three. The translations are recorded where they happen.
ITER89-P is published with density in units of 1e20 per cubic metre, which moves
a factor of `10**-0.1` into the coefficient, and with `R**1.2 a**0.3`, which
becomes `R**1.5 eps**0.3`. Both translations are stated in comments at the point
of translation, and the published form of each scaling is carried on the object
so that the translation can be checked against it.

Every exponent is verified by perturbation rather than by inspection: doubling
one input and checking that the confinement time changes by two to that exponent
goes through the same arithmetic the solver uses, which reading the dataclass
does not.

### Loss power convention as a parameter

Whether core radiation is subtracted before the confinement scaling is evaluated
is a genuine ambiguity, and it is worth a factor of 2.2 in the reported ARC gain.

The IPB98 regression was fitted with the loss power defined as total input power
less the rate of change of stored energy, with radiation not subtracted. Reactor
systems studies normally subtract core radiation and evaluate the scaling at the
power crossing the separatrix, because that is the power the edge transport
barrier actually sees, and because it is the only convention under which the
radiation terms appear explicitly in the balance at all.

Both are implemented, the separatrix convention is the default because the brief
requires the radiation terms to be visible, and the difference is reported in the
results rather than absorbed. Choosing one silently would have been the larger
error.

### Radiation, in descending order of confidence

Bremsstrahlung uses the classical `5.35e-37 Z_eff n_e**2 sqrt(T_e)` with the
relativistic and electron-electron correction of Rider (1995), which raises the
loss by 6 percent at 14 keV. The classical form is checked against a hand
evaluation and the correction can be switched off, which is what makes that check
exact rather than approximate.

Synchrotron radiation is written as the product of two factors that are handled
very differently, and the split is the point.

The optically thin emission is derived here from the Larmor formula:
`P_thin = n_e e**4 B**2 T / (3 pi eps_0 m**3 c**3)`. It is assembled from the
elementary charge, the electron mass, the electric constant, and the speed of
light rather than entered as the literal 6.2e-3 that is usually quoted, because
that literal is quoted in megawatt per cubic metre with density in units of 1e20
and is easy to transcribe with the wrong power of ten. The test suite checks the
assembled expression against the literal, which tests both at once. For ITER
parameters it comes to 1.54 MW per cubic metre, roughly fifty times the fusion
power density, which is the reason reabsorption cannot be ignored.

The escape factor is the Trubnikov asymptotic form,
`Phi = sqrt((1 - R_w) T_e B_0 / (6.04e3 a n_e20))`. The opacity coefficient in it
is the one number in the radiation model taken directly from a published
asymptotic result rather than derived here, and it is the least certain part of
the package. Isolating it in a single named constant, in a single function, with
that statement attached, is how it is handled: the uncertainty is contained
rather than distributed.

Impurity line radiation takes its loss parameter from the caller. A small table
of representative coronal-equilibrium values is provided so that the examples
have a defensible default, and its docstring says plainly that it is not a
substitute for an atomic database.

### Limits as constraints, with utilisation as the common currency

The four required limits are evaluated as numbers in the model layer and turned
into verdicts in the algorithm layer, which keeps the physics separable from the
policy about what counts as an acceptable margin.

Every check reports a utilisation, defined so that one is the boundary and above
one is a violation whichever direction the limit runs in. That is what makes an
upper bound on density comparable with a lower bound on safety factor, and it is
what allows a sweep to report which constraint binds first as a parameter is
pushed. Without a common currency, which limit binds would be a matter of
inspection rather than a computation.

A fifth constraint, the L to H transition power threshold of Martin and
colleagues (2008), was added because it is not decoration either. An H-mode
confinement scaling applied to a plasma whose loss power is below the transition
threshold is producing a number with no support. It turns out to be the binding
constraint over most of the radius sweep and at the top of every field sweep, and
it is what flags the SPARC design point as infeasible in this model. Leaving it
out would have made the high-field and large-device branches look better than the
model can justify.

### Profile factors as an optional multiplier

Fusion power is proportional to `<n**2 sigma_v(T)>`, which for any peaked profile
exceeds `<n>**2 sigma_v(<T>)`. The ratio is a pure number that depends only on
the profile shape and the mean temperature, so it can be computed once and
applied as a multiplier without giving up the zero-dimensional structure.

The default is flat, so every factor is exactly one and the model is genuinely
zero-dimensional. A parabolic family is available and its effect is reported
separately in the results. That separation is deliberate: folding a profile
factor into the default would make the model look more accurate while making the
benchmark less informative, because the reader could no longer see which part of
the discrepancy is the flat profile assumption.

The results show why a single shape cannot be adopted as a default. One shape
that brings ARC to within 12 percent of its published fusion power takes ITER to
54 percent above its own. Profile shape is a per-machine input.

### The line-averaged density, a limitation that was closed

An earlier version of these notes recorded, under "Known limitations", that the
confinement scalings are fitted against a line-averaged density while this model
supplies a volume average, and that the substitution was left uncorrected because
a correction would need the profile the model does not have. That entry has been
removed, because the correction turned out not to need a profile solution at all:
it needs only the profile shape that the caller already supplies, and for the
parabolic family it is closed form.

Both averages are integrals of the same shape function. The line average is taken
along a chord through the magnetic axis, against `d rho`; the volume average is
taken against `2 rho d rho`. Their ratio follows from the Beta function identity
`int_0^1 (1 - x**2)**a dx = sqrt(pi) Gamma(1 + a) / (2 Gamma(a + 3/2))`:

```
n_line / <n> = (1 + alpha_n) sqrt(pi) Gamma(1 + alpha_n) / (2 Gamma(alpha_n + 3/2))
```

It is exactly one for a flat density, exactly `3 pi / 8` at `alpha_n = 1/2`, and
exactly `4/3` at `alpha_n = 1`. The last two need no Gamma function and are what
the tests check the implementation against.

It is applied in the three places where a published correlation is written in
line-averaged density and nowhere else: the confinement scaling inputs, the
Greenwald limit, and the Martin L to H power threshold. At the parabolic shape
used in the results, `alpha_n = 0.4`, the ratio is 1.1446, which raises the
IPB98(y,2) confinement time by 5.69 percent at a fixed loss power, raises the
reported Greenwald fraction by the ratio itself, and raises the L to H threshold
by the ratio to the power 0.717.

What it cost. Three quantities that were previously independent of profile shape
now depend on it, which is a real loss: a reader can no longer read the Greenwald
fraction of a peaked case without knowing the shape that produced it. The peaked
ITER case now reports a Greenwald fraction of 0.959 where it reported 0.838
before, which is materially closer to a violation, and that move is a consequence
of the correction rather than of anything about ITER. The correction also assumes a chord through the
magnetic axis and nested flux surfaces of constant elongation, so a Shafranov
shifted equilibrium or a tangential chord is outside it.

What it did not cost. The default is flat, the ratio is then exactly one, and
every flat-profile number in the results is bit for bit what it was. This is a
correction that appears only when a caller has already accepted a profile shape.

What remains. The shape itself is still an input this model cannot determine, and
that is the separate "No profiles" limitation below, which stands.

### An implied confinement enhancement, to separate two failures

The fusion gain discrepancy conflates a fusion power disagreement with a
confinement disagreement, and a benchmark that reported only the gain could not
tell them apart.

Taking the published fusion power and the published auxiliary power as given, and
asking what H factor this model would need to reproduce them, separates the two
cleanly. The answer is 1.03 for ITER, 1.13 for SPARC, and 2.03 for ARC, against
assumed values of 1.0, 1.0, and 1.8. The confinement model therefore agrees with
all three sources to within 13 percent, and the gain discrepancies of 63 to 83
percent are a fusion power problem.

This diagnostic is not in the brief. It was added because without it the
benchmark reports large discrepancies and cannot say what causes them, which
makes it a table rather than a check.

## Rejected alternatives

### A one-dimensional transport solve

Solving a radial transport equation for the density and temperature profiles, and
integrating the balance over the resulting profiles, would remove the largest
single source of discrepancy in the benchmark. It would also change the character
of the repository entirely: it needs a transport model, a boundary condition at
the separatrix, a pedestal model, a time integration, and a convergence criterion,
and the answer would then depend on all of those rather than on the balance.

It was rejected on scope. The brief is a zero-dimensional benchmark and a
parametric study, and the value of a zero-dimensional model is precisely that
every number in it can be traced to an algebraic step. The cost is quantified in
the results section rather than argued about: 32 to 53 percent in fusion power,
which the parabolic profile factors move to between 12 percent low and 54 percent
high depending on the machine.

### The Albajar, Johner, and Granata synchrotron calculation

Albajar and colleagues (2001) compute the synchrotron loss from the full emission
and absorption profile rather than from a single cutoff harmonic, and their fit
is the standard reference for this term in reactor studies. It is more accurate
than the Trubnikov asymptotic form at reactor temperatures, where Trubnikov is
known to overestimate.

It was rejected on input requirements rather than on complexity. Their
formulation depends on the temperature and density profile shapes and on the
aspect ratio through a fit with many terms, and a zero-dimensional model has no
profiles to supply it. Implementing it here would mean inventing the profiles it
needs, which would make the result look more precise than its inputs.

The consequence is stated rather than hidden. The synchrotron loss at the ARC
point is 55.80 MW, nine times the bremsstrahlung loss, and it is the term that
decides whether the point produces net electricity at all. The sensitivity to the
first wall reflectivity is reported across the range 0.6 to 0.98, which spans
111.59 MW down to 24.95 MW, so a reader can see how much of the conclusion rests
on this term. If one number in this package needs replacing first, it is this
one.

### A neoclassical bootstrap calculation

The `f_BS = c sqrt(eps) beta_p` estimate is the crudest of the four required
limits. A neoclassical calculation, of the kind Sauter and colleagues give, would
compute the bootstrap current from the density and temperature gradients with the
correct collisionality dependence, and would not need a coefficient chosen by
hand.

It was rejected for the same reason as the transport solve: it needs gradients,
and a zero-dimensional model has none. The coefficient of 0.7 is the middle of
the range published estimates place it in for conventional profiles, roughly 0.5
to 1.0.

The cost of that choice is visible in the benchmark rather than concealed. The
estimate agrees with the published ARC bootstrap fraction to 2.7 percent and
overestimates ITER by 49 percent. A single coefficient applied across a factor of
seven in poloidal beta cannot do better than that, and the fact that it is right
for the high-beta case and wrong for the low-beta one is the signature of a
missing profile dependence.

### An impurity radiation database

Shipping tabulated coronal-equilibrium loss curves for a set of species would
make the line radiation term quantitative rather than indicative. It was rejected
because the constraints on this project forbid data downloads at runtime, and
transcribing a digitised curve into source would produce numbers whose provenance
a reader could not check.

Instead the loss parameter is a caller-supplied argument, a small table of
order-of-magnitude values is provided with an explicit statement of what it is
not, and the line radiation term is the smallest of the three radiation terms at
every benchmarked point.

### Fitting the profile shape per machine

The benchmark could be brought into much closer agreement by choosing a profile
shape for each machine that reproduces its published fusion power. That was
rejected because it turns a benchmark into a fit. Once the shape is chosen to
match the answer, the remaining agreement carries no information, and the
discrepancy that the reader most needs to see disappears into a free parameter.

One shape is applied to all three machines, or none is, and both cases are
reported.

### Iterating the operating point to a target gain

An alternative solver design would take a target fusion gain and solve for the
density or temperature that achieves it, which is how a design study is often
posed. It was rejected as the primary interface because it makes every result the
output of an iteration, with the reproducibility consequences described above,
and because it hides the direction of causation: in this model the gain is a
consequence of the confinement assumption, not an input to it.

A temperature solve is provided for the ignition boundary, where the question
genuinely is a root find, and it returns `None` rather than a number when its
bracket contains no sign change.

### An adaptive quadrature for the reactivity integral

The Maxwellian average of the cross-section fit could be evaluated with an
adaptive rule, which would be more accurate for a given number of function
evaluations. It was rejected because an adaptive rule stops on a floating-point
error estimate, so the number of nodes it uses depends on the arithmetic of the
machine running it, and the result is then reproducible only to the stopping
tolerance. A composite Simpson rule on a fixed grid uses the same nodes
everywhere. The same reasoning applies to the profile factor quadrature.

## Known limitations

### No profiles

There is one density and one temperature. Every volume average of a nonlinear
quantity is therefore evaluated at the mean rather than integrated, which
understates fusion power by 32 to 53 percent at the benchmarked points. The
parabolic correction factors reduce that but do not remove it, and they introduce
a profile shape that the model cannot determine for itself.

### No transport solution, and this is the dominant uncertainty

The confinement time comes from an empirical scaling, not from a transport model.
Nothing here can say why the confinement time is what it is, whether an H factor
of 1.8 is achievable, or what would change it.

This is worth restating plainly, because it is the weakest step in the chain and
everything downstream inherits it: extrapolating a confinement scaling to a
machine outside the database it was fitted to is the dominant uncertainty in any
projection of this kind. The weakness is not merely inherited, it is amplified,
by `1 / (1 - alpha_P)`, which is 3.23 for IPB98(y,2). The three scalings provided
here disagree by a factor of 7.3 in the ARC gain, and that spread is the honest
error bar on every gain and every efficiency in the results. No conclusion in this
repository about whether a device produces net electricity should be read as more
certain than that spread allows.

### No impurity transport

Impurity concentrations are inputs, uniform in space and constant in time.
Nothing models how impurities enter the plasma, how they are transported inward,
or how they accumulate on axis, which is the mechanism that determines the
tolerable high-Z concentration in a reactor. The helium ash fraction is an input
for the same reason: there is no particle balance and therefore no ash
confinement time.

### No divertor heat flux solution

The neutron power and the power crossing the separatrix are computed, and neither
is turned into a heat flux on a target. A compact high-field device has a small
major radius and therefore a short divertor wetted area, and the parallel heat
flux is one of the constraints most likely to bind first in reality. It is not
evaluated here, and no design point in this repository should be read as having
passed it.

### No magnet engineering

The toroidal field on axis is an input, with no field at the conductor, no
critical current density, no stress analysis, and no quench consideration. The
field sweeps run to 14 T on axis without asking whether a magnet could produce
it, which is the question the ARC concept exists to answer. The cryoplant appears
only as a fixed electrical load.

### No time dependence

Everything is steady state. There is no ramp, no burn control, no thermal
stability analysis of the operating point, and no pulse length. A point reported
here as satisfying every constraint has not been shown to be reachable or stable.

### Heating and current drive are not separated

The auxiliary power is one number. A steady-state device such as ARC needs
current drive whether or not it needs heating, so an ignited point in this model
draws no auxiliary power while a real one still would. The plant accounting
floors the auxiliary power at zero for that reason, which is a modelling
limitation and is marked as one at the point where the floor is applied.

### The synchrotron term is the weakest number in the model

For the reasons given above. It is the largest radiated term at every high-field
point, it decides the sign of the net electrical output at the ARC design point,
and it rests on an asymptotic formula known to overestimate. The reflectivity
sensitivity is reported so that a reader can bound the consequence.

### Fast alpha pressure is not carried

The pressure used for the beta and bootstrap calculations is thermal only. In a
burning plasma the fast alpha population contributes several percent of the total
pressure, which is part of why the computed normalised beta is 6 to 13 percent
below the published values for ARC and ITER.

### The published values are design targets, not measurements

ARC has not been built, and neither ITER nor SPARC has operated at the point
compared against here. The benchmark is a comparison between this model and three
published projections. Agreement with them is evidence that the model is
implemented consistently with the assumptions those projections used. It is not
evidence about the machines.
