# Online Appendix for *“Enforcement Policy in a Dynamic Model of Deterrence”*

This appendix provides empirical evidence, robustness testing, technical descriptions of the computational tools developed, and a guide to codebases that replicate figures in support of *”Enforcement Policy in a Dynamic Model of Deterrence”*. Links to the relevant code bases are provided where simulations are discussed.

### Notes for Online Appendix

Contents
1 Road map to the appendix 2
2 Evidence that deterrence is a dynamic process 3
2.1 Empirical signatures of a dynamic deterrence process . . . . . . . . . . . . . 3
2.2 Dynamic signatures in US and Canadian crime data . . . . . . . . . . . . . . 7
3 Robustness 12
3.1 General nature of positive feedback (CLIFFSET and DROP) . . . . . . . . . 12
3.2 Generalizing the modeling assumptions . . . . . . . . . . . . . . . . . . . . . 16
3.3 The space of history independent policies . . . . . . . . . . . . . . . . . . . . 18
3.4 The space of history dependent policies . . . . . . . . . . . . . . . . . . . . . 19
3.5 Generality of deterrence policy analysis . . . . . . . . . . . . . . . . . . . . . 24
4 Computational tools 28
4.1 Benchmarking the simulator . . . . . . . . . . . . . . . . . . . . . . . . . . . 28
4.2 Search algorithm for history dependent policies . . . . . . . . . . . . . . . . 31
5 Replicating figures and core results 34
5.1 Figures 3.1, 3.3 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 34
5.2 Figure 3.2 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 34
5.3 Figure 3.5 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 34
5.4 Figure 3.6 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35
5.5 Figure 3.7 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35
5.6 Figure 3.8 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35
5.7 Figure 3.9 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35
5.8 Figure 4.1 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36
5.9 Figure 6.1 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36
5.10 Figure 4.3 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36
6 Equilibrium selection (Pareto-dominant equilibria) 38

