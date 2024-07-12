# To-Do List

## Current Problems

### Stochastic Implementation of Compartment Variables

$S$, $E$, $P$, $I$, $A$, and $R$ are described in de Lima *et al.,* (2024) as being implemented either in a deterministic fashion (see Appendix A. Supplementary Methods - Disease Transmission Model) or in a stochastic fashion (Appendix A. Supplementary Methods - Stochastic Implementation). The functions are modelled as follows:

$$
\frac{dS}{dt} = - \lambda_t \odot S
$$

and 

$$
S_{t} = S_{t - 1} - B(S_{t - 1}, 1 - \text{exp}(\frac{-\lambda_t}{N}))
$$