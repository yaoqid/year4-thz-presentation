# Speech Transcript: Yaoqi Dong Electric Field Part

Target length: about 7 to 9 minutes for the 8-slide focused presentation.

## Slide 1: Electric Field Computation

Good morning/afternoon. My name is Yaoqi Dong, and in this presentation I will focus only on my contribution to the group project: the electric-field part of the THz emitter simulation.

The overall team project includes several mechanisms and modelling components, but my part is the link between simulated photocurrent and the far-field terahertz electric field. I will cover four connected pieces: first, the electric-field theory based on retarded potentials; second, the MATLAB implementation of the field computation; third, the angular emission validation; and finally, the custom-model results and the bipolar pulse waveform.

The main message is that the model does not only produce plots. It provides a physically consistent pipeline from current density, through retarded-time radiation, to observable THz field waveforms and angular patterns.

## Slide 2: Electric Field Theory

The electric-field calculation is based on the far-field radiation term from Jefimenko's equations. The key idea is causality. The detector does not see the current at the same instant it happens at the source. It sees the current at the retarded time, which is the observation time minus the propagation delay, \(R/c\).

In the far-field approximation, the Coulomb and induction terms fall off faster with distance, so the dominant contribution is proportional to the time derivative of the current density. This is why the equation contains \(\dot{\mathbf{J}}\), rather than just \(\mathbf{J}\).

The double cross product is also important. It projects the current acceleration onto the transverse direction, because only the field components perpendicular to the line of sight propagate as radiation.

Physically, this explains the angular behaviour of the Photo-Dember source. If the current is radially symmetric, then on-axis transverse components cancel. But when the detector is off-axis, that symmetry is broken, and a net transverse THz field remains.

## Slide 3: Electric Field Computation

This slide shows how the theory becomes a numerical solver.

The semiconductor surface is discretised into a two-dimensional source grid. Each grid cell is treated as a small radiating current element, or equivalently a Hertzian dipole contribution. The current history is stored in MATLAB as three arrays: `Jx`, `Jy`, and `Jz`, each with dimensions \(N_y \times N_x \times N_t\).

For every observation time, the solver computes the vector from each source cell to the observation point. From this it obtains \(R_{ij}\) and \(\hat{\mathbf{R}}_{ij}\), then applies the discrete version of the far-field integral.

The time derivative of current is computed using a central finite difference. This matches the implementation in `compute_THz_E_from_Jx.m`, where the derivative is calculated for all spatial points at once using MATLAB slice indexing. This is important because \(\dot{\mathbf{J}}\) directly controls the radiated waveform amplitude and shape.

The computational cost scales as \(O(N_x N_y N_t)\) for one observation point, so the implementation has to avoid unnecessary spatial loops.

## Slide 4: Retarded-Time Algorithm

The practical challenge is retarded-time interpolation.

For each source cell, the retarded time is slightly different because each cell has a different distance to the detector. That means the solver cannot simply use one common time index for the whole grid.

The algorithm has three stages. First, it differentiates the current density in time. Second, it maps the retarded time into a fractional array index. Since the retarded time normally falls between two stored time samples, the code applies linear interpolation between the two neighbouring time indices. Third, it applies the vector projection and sums the contributions from all source cells.

In the code, this is vectorised using linear indexing. The variables `idx1` and `idx2` select the correct time samples for every spatial cell. A validity mask then removes contributions whose retarded times fall outside the simulation window. This keeps the implementation both physically correct and computationally efficient.

## Slide 5: Photo-Dember Source and Angular Validation

For the Photo-Dember mechanism, the lateral current comes from the carrier-density gradient. In the standard model, the current is proportional to the spatial derivative of the Gaussian excitation profile and to the surface carrier density as a function of time.

The left plot shows a snapshot of \(|\partial \mathbf{J}/\partial t|\). The ring shape is important: the derivative is strongest at the edge of the Gaussian beam, where the carrier-density gradient is largest. At the centre, the carrier density is high, but the gradient is close to zero.

This source pattern explains the angular validation result on the right. The simulated peak-to-peak field, \(E_{pp}\), reaches its maximum at about 25 degrees. This agrees with the experimental range of roughly 20 to 30 degrees reported in the literature. The field also falls close to zero at 90 degrees after normalisation.

So the validation is not just a curve fit. The source-side current distribution and the far-field angular pattern are physically consistent with each other.

## Slide 6: Bipolar Pulse Waveform

This slide focuses on the built-in-field mechanism and the bipolar THz waveform.

In this case, the main source is a vertical drift current in the depletion field. The current is proportional to carrier density, mobility, charge, and the built-in electric field.

Because the radiated field is proportional to \(\dot{\mathbf{J}}\), a current that rises quickly and then decays more slowly naturally produces a bipolar pulse. The rising edge gives one lobe, and the decay gives the opposite lobe.

The two plots show the simulated \(E_x\) and \(E_y\) components. They both show the expected asymmetric bipolar waveform. The leading signs are opposite because the observation point is displaced off-axis in both \(x\) and \(y\), so the geometric projection gives different signs for the two transverse components.

Both components return to zero within the simulation window, which indicates that the waveform is not being truncated.

## Slide 7: Custom Built-In-Field Validation

The custom mode lets the user enter their own current-density expressions. To validate that pathway, I used a simple custom expression to reproduce the built-in-field mechanism.

In the physical built-in model, the vertical current is proportional to carrier density. Therefore, in custom mode I entered \(J_x = 0\), \(J_y = 0\), and \(J_z = -n\).

The custom expression omits physical prefactors such as charge, mobility, and field strength, so the amplitude is not expected to match. But the waveform shape should match if the parser and far-field integration are working correctly.

The two plots show that the bipolar shape and zero crossing are reproduced. This confirms that the custom expression is being converted into current-density arrays correctly, and that those arrays pass through the same far-field solver.

## Slide 8: Custom Photo-Dember Validation

Finally, I tested the custom pathway for the Photo-Dember mechanism.

Physically, the Photo-Dember current is proportional to the gradient of carrier density, with a prefactor involving charge and the difference between electron and hole diffusion coefficients. In custom mode, I entered \(J_x = dnx\), \(J_y = dny\), and \(J_z = 0\).

Again, the amplitude is not expected to match because the prefactor \(q(D_e - D_h)\) is omitted. The important comparison is the waveform shape.

The standard and custom waveforms agree in their main bipolar structure, which validates the custom parser and the far-field integrator. One difference is that the custom mode has a stronger trailing feature. This is physically meaningful: the custom mode uses the evolving carrier density \(n(x,y,t)\), whereas the standard mode factorises the current into a scalar surface density \(n_{\mathrm{surf}}(t)\) times a fixed Gaussian gradient.

So the custom mode is not just a convenience feature. It also exposes assumptions in the standard model and allows more flexible physical testing.

To conclude, my contribution was to build and validate the electric-field engine: starting from current density, applying retarded-time far-field radiation, validating angular emission, and confirming that custom current models produce physically consistent THz waveforms.
