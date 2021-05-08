# Mathematica analysis of QEdark crystal-form factors and Liquid Xenon form factors

### Datasets:
The crystal form-factors are produced by `QEdark` while the code to produce the liquid xenon form-factors is not yet public. However, a description for how to calculate the liquid xenon form-factors can be found at (https://arxiv.org/abs/1703.00910).

There are three types of datasets:
- integrated form-factors (`*.tar`): these form factors are integrated over a Maxwell-Boltzmann velocity distribution with $v_0=230$ km/s, $v_{\rm esc}=600$ km/s, $v_E=240$ km/s, and $\Delta v_E=15$ km/s. The appropriate analysis notebook is called `QEdark.nb`.
- un-integrated form-form factors (`Si_f2.dat`,`Ge_f2.dat`): to use these form factors, the user will need to specify a specific dark matter velocity profile to calculated the expected rates. The appropriate analysis notebook is called `QEdark-f2.nb`. __CAUTION__! The `*.dat` files provided were produced with an older version of `QEdark`. If you calculate crystal form-factors using the current version of `QEdark`, please remove the factor containing the k-point weighting:`wk/4`. 
- xenon form factors (`Xe_v2.mx`): these are required to calculated the dark matter-electron scattering rate in xenon. The appropriate analysis notebook is `LDM_Xenon.nb` and it requires the use of `myUnits.m`. 
