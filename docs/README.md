
# `spin-transport`

This repository contains the (developing) open-source code for simulating bulk spin transport (diffusion and separation) in solid media. Multi-spin-species simulations including magnetic resonance are in development.

## Installation

## The `spin_transport_simulation` class

### Physical constants

Constants are available via the class property `constants`, a `struct` with fields including each constant. Documentation for each should be available via the property `docs.constants`. Below, we print the documentation for each.


```matlab
field_names = fieldnames(s.constants);
for i = 1:numel(field_names)
  doc_name = ['s.docs.constants.',field_names{i}];
  print_doc(doc_name)
end
```

    .constants.ge
    	needs documentation
    .constants.gp
    	needs documentation
    .constants.hb
    	reduced Plancks constant
    	units: m^2 kg/s
    	source: https://en.wikipedia.org/wiki/Planck_constant#Value
    .constants.gamma_e
    	needs documentation
    .constants.gamma_p
    	needs documentation
    .constants.mu
    	needs documentation
    .constants.kB
    	needs documentation
    .constants.NA
    	needs documentation
    .constants.mu_B
    	needs documentation
    .constants.mu_e
    	needs documentation
    .constants.mu_N
    	needs documentation
    .constants.mu_p
    	needs documentation


### Parameters

Parameters are available via the class property `parameters`, a `struct` with fields including each parameter. Documentation for each should be available via the property `docs.parameters`. Below, we print the documentation for each.


```matlab
field_names = fieldnames(s.parameters);
for i = 1:numel(field_names)
  doc_name = ['s.docs.parameters.',field_names{i}];
  print_doc(doc_name)
end
```

    .parameters.MwPS
    	needs documentation
    .parameters.dPS
    	needs documentation
    .parameters.nAMPS
    	needs documentation
    .parameters.concDPPH
    	needs documentation
    .parameters.concPS
    	needs documentation
    .parameters.den2
    	needs documentation
    .parameters.Delta_2
    	needs documentation
    .parameters.MwDPPH
    	needs documentation
    .parameters.dDPPH
    	needs documentation
    .parameters.nAMDPPH
    	needs documentation
    .parameters.den3
    	needs documentation
    .parameters.Delta_3
    	needs documentation
    .parameters.Gamma_2
    	needs documentation
    .parameters.Gamma_3
    	needs documentation
    .parameters.grad
    	needs documentation
    .parameters.Bd_2
    	needs documentation
    .parameters.Bd_3
    	needs documentation
    .parameters.B_d
    	needs documentation
    .parameters.B0
    	needs documentation
    .parameters.B1max_p_nom
    	needs documentation
    .parameters.B1max_e_nom
    	needs documentation
    .parameters.temp
    	needs documentation
    .parameters.tPFunc
    	needs documentation
    .parameters.rPFunc
    	needs documentation
    .parameters.T12sec
    	needs documentation
    .parameters.T13sec
    	needs documentation
    .parameters.T12
    	needs documentation
    .parameters.T13
    	needs documentation
    .parameters.T22sec
    	needs documentation
    .parameters.T23sec
    	needs documentation
    .parameters.T22
    	needs documentation
    .parameters.T23
    	needs documentation
    .parameters.g
    	needs documentation
    .parameters.G
    	needs documentation
    .parameters.D
    	needs documentation
    .parameters.B_r
    	needs documentation
    .parameters.c
    	needs documentation
    .parameters.t_max_sec
    	needs documentation
    .parameters.t_max
    	needs documentation
    .parameters.r_max_nm
    	needs documentation
    .parameters.r_max
    	needs documentation
    .parameters.n_r
    	needs documentation
    .parameters.pulse
    	needs documentation
    .parameters.plot
    	needs documentation
    .parameters.rho_1_langevin
    	needs documentation
    .parameters.rho_2_langevin
    	needs documentation
    .parameters.rho_3_langevin
    	needs documentation
    .parameters.d_rho_1_langevin
    	needs documentation
    .parameters.d_rho_2_langevin
    	needs documentation
    .parameters.d_rho_3_langevin
    	needs documentation



```matlab

```
