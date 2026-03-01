This is the official repository for the paper **"Kinematic Fingerprints of Dynamical Ejection: High-Velocity Blue Straggler Stars in the Halo of NGC 188"** by Huanbin Chi.

## Overview
This repository contains the Python scripts and data processing pipelines used to analyze the dynamics of the ancient open cluster NGC 188 (~7 Gyr) using *Gaia* DR3 astrometry. 

By applying the UPMASK algorithm, we derived a highly complete membership catalog and identified a clean population of Blue Straggler Stars (BSSs). Most notably, this codebase includes the tools used to discover and trace the orbital history of two runaway BSSs in the distant cluster halo ($r > 2 r_h$), providing direct kinematic evidence for recent binary-binary scattering in the cluster core.

## Key Scientific Highlights
* **Robust Membership:** Identification of 1,652 cluster members with high completeness in the extended halo.
* **Structural Anchoring:** Calculation of the projected half-mass radius ($r_h = 5.59$ pc) and a completeness-corrected dynamical mass ($M_{tot} = 2180 M_\odot$).
* **Orbit Back-tracing:** Integration under a Plummer potential to calculate the ejection age (~1.5 Myr) and origin radius of the runaway candidates.
* **Statistical Forecasting:** Poisson-based estimation of the cluster's long-term stellar contribution to the Galactic field.
