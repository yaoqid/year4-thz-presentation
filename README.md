# Year 4 THz Electric Field Presentation

This repository contains the final presentation deck and the source material used to build it for the ELEC0118 4th Year MEng Final Project team presentation.

## Final Presentation

- `presentation/year4-thz-electric-field-presentation.pptx` - editable PowerPoint presentation generated from the project report, MATLAB code, supplied marking rubric, presentation guidance, and UCL/EEE PowerPoint template.
- `presentation/yaoqi-electric-field-part-presentation.pptx` - focused 8-slide deck for Yaoqi Dong's contribution only: Electric Field Theory, Electric Field Computation, Far-field THz Angular Emission Validation, Custom Model Results, and Bipolar Pulse Waveform.

## Repository Contents

- `input-documents/mark_scheme.pdf` - ELEC0118 final project presentation marking rubric.
- `input-documents/presentations.pdf` - presentation preparation guidance slides.
- `input-documents/presentation_template.pptx` - supplied UCL/EEE presentation template.
- `input-documents/year4_project.zip` - original project-code archive supplied for reference.
- `project-code/` - MATLAB project files used to understand and describe the simulation implementation.
- `report-source/` - LaTeX report source and figures used as the technical basis for the deck.

## Technical Story Covered By The Deck

The presentation explains an end-to-end MATLAB modelling pipeline for terahertz emission from ultrafast semiconductor photocurrents:

1. Carrier generation, diffusion, recombination, and source-current formation.
2. Retarded-potential far-field electric-field computation from `Jx`, `Jy`, and `Jz` histories.
3. Vector projection of current acceleration into transverse radiated fields.
4. Time-domain waveform results and spatial field maps.
5. Validation against the published angular emission lobe, with peak normalised amplitude around 25 degrees.
6. Critical discussion of modelling assumptions and future calibration work.

## Main Code Files Referenced

- `project-code/simulateTHzFromDiffusion.m` - carrier dynamics and current-density generation.
- `project-code/compute_THz_E_from_Jx.m` - vectorised single-observation far-field electric-field calculation.
- `project-code/compute_THz_E_spatial_2D.m` - 2D observation-grid field calculation.
- `project-code/result.m` - angular emission validation plotting data.

## Notes

The deck is designed around the presentation guidance principle of one clear proposition per slide. It also explicitly addresses the mark-scheme areas: problem and objectives, background theory, technical approach, results, critical analysis, conclusions, future work, visual organisation, and key message.
