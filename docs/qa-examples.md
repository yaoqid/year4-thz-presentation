# Q&A Examples: Yaoqi Dong Electric Field Part

## 1. Why is the radiated electric field proportional to \(\dot{\mathbf{J}}\) rather than \(\mathbf{J}\)?

A steady current does not radiate electromagnetic waves. Radiation is produced by accelerating or decelerating charge, so the far-field term depends on the time derivative of current density. In the implementation this is why `compute_THz_E_from_Jx.m` first computes `dJx`, `dJy`, and `dJz` using central finite differences before evaluating the field.

## 2. Why can you ignore the Coulomb and induction terms in Jefimenko's equations?

The detector is placed in the far-field region, where the observation distance is much larger than the source size. In that limit, the \(1/R^3\) Coulomb term and \(1/R^2\) induction term decay much faster than the \(1/R\) radiation term. The model therefore keeps the dominant radiative contribution.

## 3. What does the double cross product physically mean?

The term \(\hat{\mathbf{R}} \times [\hat{\mathbf{R}} \times \dot{\mathbf{J}}]\) projects the current acceleration onto the transverse plane. Only field components perpendicular to the propagation direction radiate to the far field. This is also why observation angle matters.

## 4. Why is the Photo-Dember emission weak on-axis?

The Photo-Dember current is radially symmetric for a symmetric excitation spot. On the optical axis, transverse contributions from opposite sides of the source cancel each other. Off-axis observation breaks this symmetry, leaving a net transverse field and producing the angular emission lobe.

## 5. Why does the angular emission peak around 25 degrees?

There is a trade-off between two effects. At small angles, the transverse projection of the radially symmetric current is weak. At larger angles, distance and projection effects reduce the observed amplitude. The simulation gives a maximum at 25 degrees, consistent with the reported experimental range of about 20 to 30 degrees.

## 6. Why is linear interpolation needed for retarded time?

The retarded time \(t_{\mathrm{ret}} = t - R/c\) usually does not fall exactly on one of the stored simulation time steps. The code maps it to a fractional array index and linearly interpolates between the two neighbouring time samples. Without this, the waveform would have time-discretisation artefacts.

## 7. Why use vectorised indexing instead of nested loops?

Each source cell has a different retarded time, so a direct implementation could require loops over \(x\), \(y\), and \(t\). That would be slow. The code constructs linear indices `idx1` and `idx2` so MATLAB can interpolate the correct retarded current derivative for all spatial cells using array operations.

## 8. What does the validity mask do?

Some retarded times fall outside the stored simulation window. The code first clamps indices to avoid array errors, then uses a validity mask to set those contributions to zero. This preserves numerical safety without adding physically invalid contributions.

## 9. Why does the built-in-field waveform become bipolar?

The built-in-field current rises quickly during carrier generation and decays more slowly during recombination. Since the radiated field follows the derivative of current, the rising phase and falling phase have opposite signs. This produces a positive-negative, or negative-positive, bipolar pulse depending on geometry.

## 10. Why do \(E_x\) and \(E_y\) have opposite leading signs in the built-in-field result?

The observation point is off-axis in both \(x\) and \(y\). The double cross-product projection maps the vertical current derivative into transverse field components. Because the observer displacement has different signs in \(x\) and \(y\), the projected \(E_x\) and \(E_y\) components can have opposite polarities.

## 11. Why does the custom built-in-field result have a different amplitude?

The custom expression \(J_z = -n\) preserves the dependence on carrier density, so it reproduces the waveform shape. It omits physical prefactors such as \(q\), \(\mu_e\), and \(E_z\), so the amplitude is expected to differ by a scale factor.

## 12. Why does the custom Photo-Dember mode show a stronger trailing feature?

The standard Photo-Dember model factorises the lateral current into \(n_{\mathrm{surf}}(t)\) times a fixed spatial Gaussian gradient. The custom expression uses the evolving carrier-density gradient \(dnx\) and \(dny\), so the spatial current pattern changes as carriers diffuse. That extra time variation can enhance the trailing part of the radiated waveform.

## 13. Does the custom mode make the standard model wrong?

No. The standard model is a simplified and calibrated pathway that is computationally efficient. The custom model is useful because it exposes the effect of relaxing some assumptions, such as using an evolving spatial gradient rather than a fixed Gaussian gradient.

## 14. What is the main limitation of your validation?

The angular validation compares normalised peak-to-peak amplitude rather than absolute calibrated THz power. It validates the shape and peak angle of the emission pattern, but absolute amplitude would require detector response, collection optics, material parameters, and calibration to be included more carefully.

## 15. What would you improve next?

I would add quantitative calibration against measured time-domain waveforms, include detector response and refraction effects, and run sensitivity studies for beam waist, carrier lifetime, mobility, depletion-field strength, and observation geometry.

## 16. How do you know the equations in the slides match the MATLAB code?

The far-field equation maps directly to `compute_THz_E_from_Jx.m`: the code computes \(R\), \(\hat{\mathbf{R}}\), \(\dot{\mathbf{J}}\), retarded-time interpolation, the transverse projection, division by \(R\), and multiplication by \(\Delta x \Delta y\). The current-source equations map to `simulateTHzFromDiffusion.m`, where `Jx_all`, `Jy_all`, and `Jz_all` are assembled before calling the far-field solver.

## 17. Why is the source grid stored as \(N_y \times N_x \times N_t\), not \(N_x \times N_y \times N_t\)?

This follows MATLAB matrix convention. `meshgrid(x,y)` returns arrays whose rows correspond to \(y\) and columns correspond to \(x\), so the current histories are stored as \(N_y \times N_x \times N_t\). The linear indexing in the field solver is written consistently with that ordering.

## 18. If the amplitude is not calibrated in custom mode, why is custom mode useful?

It is useful for testing physical dependencies and waveform shapes. If a custom current has the correct spatial and temporal dependence, the far-field solver should produce the expected waveform shape. This makes custom mode a validation and exploration tool even when absolute amplitude needs separate prefactors.
