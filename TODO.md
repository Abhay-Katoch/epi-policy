# To-Do List

## Current Problems

### Stochastic Implementation of Compartment Variables

$S$, $E$, $P$, $I$, $A$, and $R$ are described in de Lima *et al.,* (2024) as being implemented either in a deterministic fashion (see Appendix A. Supplementary Methods - Disease Transmission Model) or in a stochastic fashion (Appendix A. Supplementary Methods - Stochastic Implementation). The functions are modelled as follows:

$$
\frac{dS}{dt} = - \lambda_t \odot S
$$

and 

$$
\boldsymbol{S_{t}} = \boldsymbol{S_{t - 1}} - B(\boldsymbol{S_{t - 1}}, 1 - \text{exp}(\frac{-\boldsymbol{\lambda_{t - 1}}}{\boldsymbol{N}}))
$$