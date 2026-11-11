# Polishing pass — open items (first draft complete 2026-11-11)

State: all 56 notebooks written, executed without errors, QA clean; book builds with 0 warnings / 0 dead links.
Build commands: see RESUME_NOTES.md.

## A. TODO (polish) markers in the notebooks
| Notebook | Item |
|---|---|
| P1.5 | EPI T2* blurring / echo-train PSF demo; g-factor map; susceptibility-distortion + Gibbs demo; fact-check IEC numbers |
| B3 | passive-shim simulation only reaches ~27 ppm pk-pk (real magnets: few ppm) — finer grid / LP fit |
| B5 | birdcage mode spectra + ladder schematic; coupled-loop resonance splitting + preamp decoupling + Roemer combination; fact-check absolute SNR estimate; 7 T toy shows darker centre (should be central brightening) + RF-shimming demo |
| B6 | verify IEC rheobase/chronaxie + cardiac limit; calibrate SAR cylinder model, add local SAR; fact-check limits table, add spatial-field-gradient map |
| C1 | widget explorer + spin-echo diagram; text says min TE≈4 ms/TR≈9 ms, code computes 3.66/9.38 ms |
| C2 | 2D encoding demo + widget explorer |
| D4 | confirm sequences of Pauly & Gold / Yarnykh example images; 2nd phantom from contrast_phantom.m; VFA T1 at realistic noise; B1+ error demo |
| E3 | proper band-edge-based ripple metric |
| H3 | learned unrolled toy (24.5 % vs 5.8 % error) — improve |
| I2 | opposed-phase dark line too faint |
| J1 | slew-limited spiral (gradient jumps from 0), EPI-trajectory example, disk target |
| J2 | flyback design; EPI slice with/without SPSP showing fat shift |
| J3 | RF-along-y snapshots; time slider |
| J4 | bSSFP efficiency; Monte-Carlo g-factor check; B0-scaling example |
| J5 | compact SPIRiT implementation |

## B. Size / presentation
- Large notebooks (figure dpi / image count): I2 6.8 MB, H3 5.7 MB, F3 5.6 MB, F2 3.3 MB, I1 3.2 MB, RM6 3.2 MB.
- Clipped titles in I1/I2 figures.
- Figures never visually checked: B3, B5, J1 captions, J4, D4, later D2 sections, RP4.
- Widget ("Try it yourself") cells are not executed at build time (headless kernels hang on them); they show a note instead.

## C. Facts to verify (flagged by writers as from memory / approximate)
- Safety: IEC 60601-2-33 SAR, dB/dt (20 T/s, 0.36 ms), acoustic limits (B6, P1.5, Glossary).
- Hardware specs: vendor gradient table (B4), helium volumes of sealed magnets, "~1 M$/T", Iseult stored energy, magnet currents (B3), tissue ε_r/σ (B5, RP2, B6).
- Relaxation/physiology: 3 T T1/T2 table (D1/D2; original lists WM T2 > GM T2), BOLD Δχ and ΔR2*, venous T2*, relaxivities (D2, A3), six-peak fat model (A3), susceptibility table (A3).
- Vendor sequence-name tables (H2, Glossary) — esp. Canon names, GE fat-sat names, Philips SmartSpeed.
- RP2 order-of-magnitude estimates (33 nV/mm³ signal, 60 nV noise, 4 V/m, 2 W/kg), G2 field-strength constants.
- Literature citations written from memory throughout (check each "References" section).
- Software URLs in Software_Resources (twixtools, MRzero-Core, JEMRIS).

## D. Consistency
- Coil sensitivity convention C ∝ Bx − iBy (F3, H3, J4, J5) vs. B1⁻ sign conventions — state once in Notation and link.
- RM1 uses τ for the time constant (master notation: τ = integration variable).
- G1 `cartesian_design` vs G3 `protocol` calculators overlap — unify (maybe move to mri_tools).
- Notation_and_Symbols "first appears in" links were guessed before chapters existed — verify.
- Candidate mri_tools additions: wire_field, loop_field_exact, trapezoid / k_from_gradients, dipole_field, spinor (Cayley–Klein) simulator, dzrf.
- mt.flip_angle now uses the signed integral (fixed after A2/D2 were first written; both rebuilt).
