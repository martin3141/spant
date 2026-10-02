# Basis simulation

## Basis simulation

Basis simulation is necessary step for modern MRS analysis and the this
vignette will explain how to achieve this with spant. It is advisable to
follow the examples given in the [metabolite simulation
vignette](https://martin3141.github.io/spant/articles/spant-metabolite-simulation.md)
before following this guide.

Load the spant package:

[`library`](https://rdrr.io/r/base/library.html)`(`[`spant`](https://spantdoc.wilsonlab.co.uk/)`)`

A basis set is a collection of signals to be fit to the MRS data. In
spant we start with a list of molecular definitions containing the
relevant information for each signal - such as chemical shifts and
j-coupling values:

`mol_list`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`[`get_mol_paras`](https://martin3141.github.io/spant/reference/get_mol_paras.md)`(``"lac"``)``,`` `` `[`get_mol_paras`](https://martin3141.github.io/spant/reference/get_mol_paras.md)`(``"naa"``)``,`` `` `[`get_mol_paras`](https://martin3141.github.io/spant/reference/get_mol_paras.md)`(``"cr"``)``,`` `` `[`get_mol_paras`](https://martin3141.github.io/spant/reference/get_mol_paras.md)`(``"gpc"``)``)`

In the next step we convert these chemical properties into a collection
of signals (a spant `basis_set` object) with the `sim_basis` function.
When fitting, the signal parameters (e.g. sampling frequency) and pulse
sequence (e.g. echo-time) must match the MRS data acquisition protocol.

`basis`` ``<-`` `[`sim_basis`](https://martin3141.github.io/spant/reference/sim_basis.md)`(``mol_list``, pul_seq ``=`` ``seq_slaser_ideal``,`` `` acq_paras ``=`` `[`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)`(``N ``=`` ``2048``, fs ``=`` ``2000``, ft ``=`` ``127.8e6``)``,`` `` TE1 ``=`` ``0.008``, TE2 ``=`` ``0.011``, TE3 ``=`` ``0.009``)`` `` `[`stackplot`](https://martin3141.github.io/spant/reference/stackplot.md)`(``basis``, xlim ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``4``, ``0.5``)``, y_offset ``=`` ``50``, labels ``=`` ``basis``$``names``)`

![](spant-basis-simulation_files/figure-html/unnamed-chunk-4-1.png)

In 1H MRS broad resonances from lipids and macromolecules are often
included in addition to metabolites:

`mol_list_mm`` ``<-`` `[`append`](https://rdrr.io/r/base/append.html)`(``mol_list``, `[`list`](https://rdrr.io/r/base/list.html)`(`[`get_mol_paras`](https://martin3141.github.io/spant/reference/get_mol_paras.md)`(``"MM09"``, ft ``=`` ``127.8e6``)``)``)`` `` ``basis_mm`` ``<-`` `[`sim_basis`](https://martin3141.github.io/spant/reference/sim_basis.md)`(``mol_list_mm``, pul_seq ``=`` ``seq_slaser_ideal``,`` `` acq_paras ``=`` `[`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)`(``N ``=`` ``2048``, fs ``=`` ``2000``, ft ``=`` ``127.8e6``)``,`` `` TE1 ``=`` ``0.008``, TE2 ``=`` ``0.011``, TE3 ``=`` ``0.009``)`` `` `[`stackplot`](https://martin3141.github.io/spant/reference/stackplot.md)`(``basis_mm``, xlim ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``4``, ``0.5``)``, y_offset ``=`` ``50``, labels ``=`` ``basis_mm``$``names``)`

![](spant-basis-simulation_files/figure-html/unnamed-chunk-5-1.png)

Note the field strength is often required to simulate these broad
resonances as their linewidth is usually specified in ppm. spant also
includes the functions `sim_basis_1h_brain` and
`sim_basis_1h_brain_press` to produce commonly used sets of basis
signals:

`basis`` ``<-`` `[`sim_basis_1h_brain`](https://martin3141.github.io/spant/reference/sim_basis_1h_brain.md)`(``)`` `[`stackplot`](https://martin3141.github.io/spant/reference/stackplot.md)`(``basis``, xlim ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``4``, ``0.5``)``, y_offset ``=`` ``20``, labels ``=`` ``basis``$``names``)`

![](spant-basis-simulation_files/figure-html/unnamed-chunk-6-1.png)

Basis sets can be exported for use with LCModel with the `write_basis`
function, and sim_basis_1h_brain has the option `lcm_compat` to remove
signals that are usually generated within the LCModel package:

`lcm_basis`` ``<-`` `[`sim_basis_1h_brain`](https://martin3141.github.io/spant/reference/sim_basis_1h_brain.md)`(``)`` `[`stackplot`](https://martin3141.github.io/spant/reference/stackplot.md)`(``lcm_basis``, xlim ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``4``, ``0.5``)``, y_offset ``=`` ``20``, labels ``=`` ``basis``$``names``)`

![](spant-basis-simulation_files/figure-html/unnamed-chunk-7-1.png)
