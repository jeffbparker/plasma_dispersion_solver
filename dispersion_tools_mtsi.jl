"""
Solving plasma dispersion relations and determining stability.  For:
    --Modified two-stream instability (MSTI)

Calculating zeros of a dispersion function D(omega) or dielectric function epsilon(omega)
* Solve for local zeros using an iterative root finder (e.g., fixed-point or Newton-like)
* For analytic D, find zeros within a certain contour with the Davies method
"""

module dispersion_tools_mtsi
import SpecialFunctions
using FromFile
@from "PlasmaDispersionFunction.jl" import PlasmaDispersionFunction

# convenient aliases for functions related to plasma dispersion function Z(z)
Z = PlasmaDispersionFunction.Z
Zprime = PlasmaDispersionFunction.Zprime
Zprimehat = PlasmaDispersionFunction.Zprimehat

function D_MTSI(omegahat, kxhat, kzhat, mu, tau, ui, sigma)
    """MTSI: Kinetic dispersion relation, electrostatic

    Dimensionless variables:
    omegahat = omega / omega_LH
    kxhat= kx * U / omega_LH
    kzhat = kz * U / omega_LH
    mu = m_i / m_e
    tau = T_e / T_i
    ui = U / w_Ti
    sigma = Omega_e^2 / omega_pe^2
    """
    khat = sqrt(kxhat^2 + kzhat^2)
    ue = ui / sqrt(mu*tau)
    vteSq_over_USq = tau*mu / (2*ui^2)
    omegaLHSq_over_OmegaeSq = 1 / (mu * (sigma+1))
    omegaLHSq_over_omegapeSq = sigma / (mu * (sigma+1))
    kSq_times_lambdaDeSq = khat^2 * vteSq_over_USq * omegaLHSq_over_omegapeSq
    
    zeta_e = omegahat / (kzhat / ue)
    zeta_i = (omegahat - kxhat) / (khat / ui)
    lambda = kxhat^2 * vteSq_over_USq * omegaLHSq_over_OmegaeSq

    # perform evaluation of special functions.  Note that Zhat(z) = -Z''(z) / 2
    Ixval = SpecialFunctions.besselix(0, lambda)   # besselix(nu, x) = exp(-|Real(x)|) * I_nu(x); prevents overflow at large x
    Zval_e, Zprimeval_e, Zhatval_e = Zprimehat(zeta_e)
    Zval_i, Zprimeval_i, Zhatval_i = Zprimehat(zeta_i)

    # compute terms going into dispersion function
    Poisson_term = kSq_times_lambdaDeSq   # intrinsic Poisson response.  Term is neglected to assume quasineutrality
    chi_i = -tau/2 * Zprimeval_i
    chiprime_i = tau / (khat / ui) * Zhatval_i

    chi_e1 = 1 - Ixval  # precision is okay as long as lambda > 1e-10
    chi_e2 = -Ixval * Zprimeval_e / 2
    chiprime_e2 = Ixval / (kzhat / ue) * Zhatval_e
    D = Poisson_term + chi_i + chi_e1 + chi_e2
    Dprime = chiprime_i + chiprime_e2
    return D, Dprime
end

end