# Scientific status and implementation plan

## Implementation plan

1. **Completed:** replace the article layout with a thesis-style `book` project containing a
   title page, abstract, notation, eight numbered chapters, appendices, and a
   reproducibility record.
2. **Completed at analytical-framework level:** expand the derivations from linear elasticity through Rayleigh--Lamb
   eigenfields, physical thermoelastic forcing, causal modal response, and
   moving/periodic sources. Mark unresolved normalization or experimental inputs
   explicitly instead of hiding them behind unit overlap.
3. **Completed for verified dispersion and oscillator checks:** refactor the preliminary monolithic numerical script into importable modules
   for materials, dispersion, causal response, plotting, and validation, with one
   reproducible entry point and automated tests.
4. **Completed:** preserve the existing scalar-overlap outputs under a clearly labeled
   `results/preliminary/` area. They are framework checks only, not physical mode
   purity or absolute thermoelastic predictions.
5. **Completed for available tests:** verify analytic limits, characteristic residuals, root-grid convergence, and
   zero initial conditions; compile the complete thesis and check references.
6. **Pending final publication step:** publish the result on a new branch with a draft pull request that reports
   scientific limitations as prominently as numerical achievements.

## Verified results

- The isotropic Rayleigh--Lamb characteristic equations are solved numerically
  for representative aluminum properties recorded in the calculation metadata.
- The current verified run has maximum scaled determinant residual
  `2.64e-10`, scan-density difference `4.11e-13`, A0 low-frequency median error
  `0.136%`, and S0 error `0.0105%`.
- The first computed symmetric higher-order minimum is near
  `fd = 2.8568 MHz mm`; it is retained as a numerical benchmark only after a
  primary-source comparison is recorded in the thesis.
- The low-frequency numerical branches converge to the plane-stress extensional
  `S0` limit and Kirchhoff--Love flexural `A0` limit.
- The forced-oscillator implementation satisfies zero displacement and velocity
  at switch-on to floating-point precision.

## Preliminary results

- Existing settling-time, `N95`, and time-dependent two-mode ratios use a scalar
  unit-overlap source model and a prescribed constant damping ratio.
- The earlier frequency window is an illustrative signal-processing choice, not
  a detector model. It must not be tuned or interpreted as a physical prediction.
- Existing high-order root labels are sorted-root labels; continuous eigenfield-
  based branch tracking remains to be completed before high-order dispersion is
  treated as definitive.

## Unresolved derivations

- Full normalized Lamb eigenfields and a numerically stable modal-mass convention
  across propagating, cutoff, and ZGV regimes.
- Equivalence of body-force, surface-traction, and weak-work forms for the chosen
  thermoelastic boundary model.
- Uniform asymptotics for coalescing stationary points and endpoints.
- A detector-aware settling functional and a burst metric tied to a declared
  observable (amplitude, energy, or demodulated envelope).

## Missing physical inputs from the author/experiment

- Aluminum alloy and measured elastic constants, density, thermal properties,
  optical reflectivity, absorption depth, and thermal boundary conditions.
- Plate thickness tolerance and lateral geometry.
- Laser wavelength, absorbed power/fluence, beam convention and width, temporal
  waveform, scan trajectory, and damage-regime confirmation.
- Frequency-dependent attenuation and detector transfer function/bandwidth.

## Proposed next calculations

1. Implement eigenvector continuation using normalized displacement/stress
   overlap between adjacent wavenumber points.
2. Compute thermoelastic weak-work overlaps for finite absorption depth and test
   symmetric/antisymmetric selection rules.
3. Replace the toy scalar response with those overlaps and a declared detector
   filter, then recompute settling time and burst length.
4. Treat ordinary, cutoff, and ZGV regimes separately and compare numerical
   transients with the appropriate asymptotic forms.
