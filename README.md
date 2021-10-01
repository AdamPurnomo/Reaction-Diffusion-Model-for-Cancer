# Reaction-Diffusion-Model-for-Cancer

### Overview

A Python implementation of reaction-diffusion model simulation for invasive cancer. The model is based on this [study](https://cancerres.aacrjournals.org/content/canres/56/24/5745.full.pdf) <sup>[1]</sup> where the cancer growth is facilitated by the diffusion of excess acid production from cancer cells to the surrounding healthy tissue and described with reaction and diffusion equation. Explicit finite difference method is used to solve the differential equation and simulate the cancer growth.

#### Reference

<sup>[1]</sup> Gatenby RA, Gawlinski ET. (1996). A reaction-diffusion model of cancer invasion.
Cancer Research, 56(24):5745-53. PMID: 8971186.

### Examples

There are 3 treatment scenarios simulated in this project: with chemotherapy, with surgery, and with a combination of chemotherapy and surgery.

#### Surgery

<p align="center">
  <img width=50% height=50% src="https://github.com/AdamPurnomo/Reaction-Diffusion-Model-for-Cancer/blob/main/With%20Surgery.gif?raw=true">
</p>

#### Cemotherapy 

<p align="center">
  <img width=50% height=50% src="https://github.com/AdamPurnomo/Reaction-Diffusion-Model-for-Cancer/blob/main/With%20Chemo.gif?raw=true">
</p>

#### Surgery and Chemotherapy 

<p align="center">
  <img width=50% height=50% src="https://github.com/AdamPurnomo/Reaction-Diffusion-Model-for-Cancer/blob/main/With%20Surgery%20and%20Chemo.gif">
</p>

### Dependencies

* numpy 1.19.2
* matplotlib 3.3.4




