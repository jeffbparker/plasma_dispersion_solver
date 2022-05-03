# Solving kinetic plasma dispersion relations in Julia
This package provides some libraries and notebooks to aid in solving kinetic plasma dispersion relations in Julia.  We solve the Vlasov-Poisson system, with cases including either magnetized and unmagnetized electrons.

## Packages
* Package `PlasmaDispersionFunction` contains some functions for evaluating the plasma dispersion function $Z(z)$ as well as its derivatves $Z'$ and $Z''$.  These are robust in the sense that arguments that would otherwise numerically be troublesome are handled specially using the asymptotic series.  In particular, the exact expression $Z'(z) = -2[1 + zZ(z)]$ suffers truncation error at large $z$, and similarly for $Z''$.  
* Package `continuation` implements some helpers for doing numerical continuation in a single parameter.  This is useful when repeatedly solving a dispersion relation as a parameter varies over some range, and the solution is expected to be continuous along that branch.
* Package `plasma_dispersion_tools` and `dispersion_tools_mtsi` provide functionality for computing the dispersion relation residual $D(\omega)$ for various instabilities or assumptions, where $D(\w)=0$ indicates a frequency $\omega$ that corresponds to an eigenmode.
* Package `tools` implements an algorithm from complex analysis using the argument principle that enables one to count and find the zeros of complex analytic functions.  This package implements ideas contained in:
	* Delves & Lynesse, A Numerical Method for Locating the Zeros of an Analytic Function, Math. Comp., 21, 543 (1967).
	* B. Davies, Locating the Zeros of an Analytic Function, J. Comp. Physics, 66, 36-49 (1986).
	* The squircle contour is discussed in C. D. Stephens et al., Quasilinear gyrokinetic theory: a derivation of QuaLiKiz, J. Plasma Phys. (2021) 87, 905870409
* There are some corresponding packages prefixed by `big_` which performs the dispersion calculations in arbitrary precision arithmetic using a combination of `BigFloat` and the `ArbNumerics` package.  (These `big_` packages are somewhat incomplete)

## Notebooks
In this repository are several notebooks that use the packages to actually carry out the task of finding the zeros of the dispersion relation and therefore finding the eigenmodes.  Solved are the Ion-Ion Streaming Instability and several variations of the Modified Two-Stream Instability.
