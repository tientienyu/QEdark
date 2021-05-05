# QEdark

## The `QEdark` code is based on work as published in [arXiv:1509.01598](https://arxiv.org/abs/1509.01598) by Rouven Essig, Marivi Fernandez-Serra, Jeremy Mardon, Adrian Soto, Tomer Volansky, and Tien-Tien Yu

The following repository contains the analysis notebooks for [QEdark](https://github.com/adrian-soto/QEdark_repo), which are needed for calculating dark matter-electron scattering in semiconductor targets. We provide both Mathematica and Jupyter notebooks, which utilize the tabulated form-factors calculated using [QEdark](https://github.com/adrian-soto/QEdark_repo). 

### Datasets:
There are three types of datasets:
- integrated form-factors (silicon, germanium): these form factors are integrated over a Maxwell-Boltzmann velocity distribution with $v_0=230$ km/s, $v_{\rm esc}=600$ km/s, $v_E=240$ km/s, and $\Delta v_E=15$ km/s. The appropriate analysis notebook is called `QEdark.nb`.
- un-integrated form-form factors (silicon, germanium): to use these form factors, the user will need to specify a specific dark matter velocity profile to calculated the expected rates. The appropriate analysis notebook is called `QEdark-f2.nb`. 
- xenon form factors: these are required to calculated the dark matter-electron scattering rate in xenon. The appropriate analysis notebook is `LDM_Xenon.nb` and it requires the use of `myUnits.m`. 
