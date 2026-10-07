# tools/aao_rc_point — radiative correction at a fixed kinematic point

Computes RC = sigma_obs/sigma_Born of aao_rad's Mo-Tsai machinery at a
FIXED point (Q2, xB or W, t, phi): no event generation, no binning, no
migration. Built for point-by-point comparison with EXCLURAD delta(phi).

## Usage
```
./build_rc_point.sh
./scan_phi.py -E0 5.75 -Q2 1.125 -xb 0.1372 -t 0.12 -vcut 0.2
```
Output: rc_point_scan.txt (phi, RC). ~17 s per phi point at 2e7 tries.

## How it works
aao_rc_point.F is a surgical copy of aao_rad.F (regenerate with
make_rc_point.py after any physics change in aao_rad.F!) with:
- electron kinematics frozen via degenerate Q2/E' input ranges;
- hadronic CM angles (cos theta*, phi) frozen from point.inp;
- on every MC try TWO weights accumulate through the identical code
  path: the observed one (with radiative factors) and a Born one
  (soft branch, radiative factors set to 1; the soft sampling measure
  integrates to unity). RC = sum(obs)/sum(born): all jacobians,
  region factors and common constants cancel exactly by construction;
- event building/output is bypassed; fixed number of tries.

## Caveats
- Hadron angles are fixed at the VERTEX; EXCLURAD fixes the OBSERVED
  proton t/phi. Within the vcut window the difference is small but
  nonzero for hard photons.
- aao_rad's Mo-Tsai RC gives a FLATTER phi-dependence than the exact
  Bardin-Shumeiko integral (EXCLURAD), and sits below it.  How much is
  energy dependent, so quote the point: at the 5.75 GeV CLAS6 pi0 test
  point (2026-07) the phi-dependence came out ~10x flatter; at Q2 = 1.5,
  xB = 0.30, -t = 0.30, v < 0.2, E0 = 10.6 GeV (2026-10) it is only 1.6x
  flatter -- swing 0.061 against 0.097 -- while the offset is 10%.
- That offset is not the structure functions: pi0.vpk2021, pi0.amp2609 and
  pi0.vpk2026 agree with each other to 0.5%.  One piece of it is identified
  and is a genuine omission -- aao_rad implements no vacuum polarization at
  all (no vacuum/vpol/delvac anywhere in aao_rad/*.F), while EXCLURAD
  carries delta_vac = +3.6%, phi independent.  The remaining ~6% is NOT
  decomposed.  It is not the peaking approximation: aao_rad integrates the
  photon solid angle exactly, by importance sampling in 5 regions
  (aao_rad.F:215 ff) -- two narrow windows on the incident and scattered
  electron directions, two shells, and a fifth region covering the full
  cos(theta_k) and phi_k with the peak windows subtracted, each carrying its
  compensating mcfac/mpfac weight.  Open candidates: the soft/hard split and
  what is exponentiated, the vertex-versus-observed hadron angles (above),
  and the MAID/DVMP seam at W = 1.8 inside the radiative integrand.
- Until that is settled, use EXCLURAD for RC numbers; use aao_rad for event
  samples (acceptance, MM2 shapes).


# tools/check_lt_convention — are the two Born treatments consistent?

`dvmpw` returns both a cross section `sigma0` and the structure functions.
`aao_norad` samples events from `sigma0`; the radiative branch of `aao_rad`
(`aao_rad.F:977`) instead rebuilds the Born from the structure functions in
the AO convention.  The two must agree, and this program checks that they
do at a fixed kinematic point.

```
gfortran -fno-automatic -ffixed-line-length-none \
  tools/check_lt_convention.F aao_rad/dvmpw.F aao_rad/dvmpx.F \
  -o tools/check_lt_convention
tools/check_lt_convention                      # defaults: 10.6 1.5 0.30 0.30
tools/check_lt_convention 5.75 1.125 0.1372 0.12   # E0 Q2 xB |t|
```

Prints `ratio = born_AO/sigma0` over phi and exits non-zero if it strays
from 1 by more than 1e-4.  Before `c83ce8e` it ran 0.801 / 1.000 / 1.987 at
phi = 0 / 90 / 180 -- exactly 1 at 90 and 270, where `cos(phi) = 0` removes
the LT term, which is the signature of a convention mismatch in `sigma_LT`
rather than anything radiative.
