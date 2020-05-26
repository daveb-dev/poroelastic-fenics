# poroelastic-fenics

This repo is about the development of single-phase poroelastic models using a mixed finite element. The mass and momentum conservation equations are discretized into weak form considering flow velocity as an extra primary variable of the problem. The discretized weak form is implemented in FeniCS.

# Validation of the model:
- 1D Tarzhagi's consoldiation problem - https://github.com/rksin8/poroelastic-fenics/tree/master/terzaghi
- 2D Mandel's consolidation problem - https://github.com/rksin8/poroelastic-fenics/tree/master/mandel

# Real case example
- Simuluation of heterogeous reservoir SPE10 data - https://github.com/rksin8/poroelastic-fenics/tree/master/example03_hetero_iso
- Simulation of layered reservoir - https://github.com/rksin8/poroelastic-fenics/tree/master/example01
