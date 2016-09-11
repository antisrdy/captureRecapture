# Capture Recapture

Assessment realized in 2016 on capture recapture.

Definition of caputre recapture from Wikipedia: *"A portion of the population is captured, marked, and released. Later, another portion is captured and the number of marked individuals within the sample is counted. Since the number of marked individuals within the second sample should be proportional to the number of marked individuals in the whole population, an estimate of the total population size can be obtained by dividing the number of marked individuals by the proportion of marked individuals in the second sample."*

Considering *p* as the probability of being captured, *M* as the population size, and *x* as available data, the present work aims at:
* Computing theoritical distributions of conditionnal probabilities of *p/M,x* and *M/p,x*
* Find a way to estimate them (Gibbs-Sampler)
* Estimating Monte-Carlo error

Theoritical answers are inspired by (https://www.ceremade.dauphine.fr/~xian/bcs/bcap.pdf)
