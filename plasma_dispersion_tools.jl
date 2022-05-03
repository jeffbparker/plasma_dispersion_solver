"""
Solving plasma dispersion relations and determining stability.  For:
    --ion-ion streaming instability (Foote and Kulsrud)
    --alpha ion-ion streaming instability

Calculating zeros of a dispersion function D(omega) or dielectric function epsilon(omega)
* Solve for local zeros using an iterative root finder (e.g., fixed-point or Newton-like)
* For analytic D, find zeros within a certain contour with the Davies method

Other methods of determining stability
* Count zeros of analytic epsilon within a contour using the argument principle (or use real line + small imag part)
    --This is essentially the Nyquist diagram
* Use the Penrose criterion to determine electrostatic instability.
    --This method does not depend on the wavenumber, AND one can easily substitute non-Maxwellian distribution functions
"""

module plasma_dispersion_tools
import SpecialFunctions

function Z(z)
    """Plasma dispersion function Z(z).
    Faddeeva function is w(z) = erfcx(-iz).
    Plasma dispersion function is Z(z) = i * sqrt(pi) * w(z)
    """
    out = im * sqrt(pi) * SpecialFunctions.erfcx(-im*z)
    return out
end

function Zprime(z)
    """Derivatve Z'(z).  Note that at large |z| there can be precision loss, so multiple precision helps"""
    out = -2 * (1 + z * Z(z))
    return out
end

function Zhat(z)
    """Derivative d(zZ)/dz = zZ' + Z"""
    Zval = Z(z)
    Zprimeval = -2 * (1 + z * Zval)
    out = z * Zprimeval + Zval
    return out
end

function Zprimehat(z)
    """Return Z(z), Zprime(z), and Zhat(z) all together."""
    Zval = Z(z)
    Zprimeval = -2 * (1 + z * Zval)
    Zhatval = z * Zprimeval + Zval
    return Zval, Zprimeval, Zhatval
end

function D_foote_real(omegahat_i, kparallelhat, kperphat, ui, tau)
    """Use dispersion function in Foote & Kulsrud for ion-ion streaming instability.
        
        --D(omegahat) in Foote neglects the intrinsic Poisson response (i.e. assumes quasineutrality)
        --Also scales the dielectric function epsilon by a factor to yield D.
        --Here: assume omega is imaginary: and give the imaginary part (a real quantity)
        --Then D is real, so return a real number.
    
    Here, evaluate both D(omegahat_i) for adiabatic electrons, Maxwellian ions.

    
    Dimensionless variables:
    omegahat_i = imag(omega) / omega_pi
    kparallelhat = k_parallel * lambda_De
    kperphat = k_perp * lambda_De
    ui = U / w_Ti
    tau = T_e / T_i
    """
    tau2 = 2 / sqrt(tau)
    khat = sqrt(kparallelhat^2 + kperphat^2)

    zeta_iplus = (im*omegahat_i + kparallelhat * ui * tau2) / (khat * tau2)
    Zval_iplus = Z(zeta_iplus)

    D = 1/tau + 1 + zeta_iplus * Zval_iplus
    Dr = real(D)
    return Dr
end

function D_aiisi_real(omegahat_i, kparallelhat, kperphat, nu_d, nu_alpha, tau_e, tau_alpha, ui)
    """dispersion function for alpha ion-ion streaming instability.
        
        --D(omegahat) neglects the intrinsic Poisson response (i.e. assumes quasineutrality)
        --Assume adiabatic electron response
        --Also scales the dielectric function epsilon by a factor to yield D.
        --Here: assume omega is imaginary: and give the imaginary part (a real quantity)
            -->This last assumption means that D is real, so return a real number.

    Dimensionless variables:
    omegahat_i = imag(omega) / omega_pd
    kparallelhat = k_parallel * lambda_De
    kperphat = k_perp * lambda_De
    nu_d = n_d / n_e
    nu_alpha = n_alpha / n_e
    tau_e = T_e / T_d
    tau_alpha = T_alpha / T_d
    ui = V / w_Td
    """
    tau_e2 = sqrt(2 / (nu_d * tau_e))
    tau_alpha2 = sqrt(tau_alpha / (nu_d * tau_e))
    khat = sqrt(kparallelhat^2 + kperphat^2)

    zeta_d = im*omegahat_i / (khat * tau_e2)
    zeta_t = im*omegahat_i / (khat * sqrt(2/3) * tau_e2)
    zeta_aplus = (im*omegahat_i + kparallelhat * ui * tau_e2) / (khat * tau_alpha2)

    Zval_d = Z(zeta_d)
    Zval_t = Z(zeta_t)
    Zval_aplus = Z(zeta_aplus)

    # scaled susceptibilities chi --- multiplied by k^2 lambda_Dd^2
    chi_e = 1 / (nu_d * tau_e)
    chi_d = 1 + real(zeta_d * Zval_d)
    chi_t = 1 + real(zeta_t * Zval_t)
    chi_alpha = 4 * nu_alpha / (nu_d * tau_alpha) * 2 * real(1 + zeta_aplus * Zval_aplus)
    Dr = chi_e + chi_d + chi_t + chi_alpha
    return Dr
end


function D_aiisi_electrons_real(omegahat_i, kparallelhat, kperphat, nu_d, nu_alpha, tau_e, tau_alpha, ui)
    """dispersion function for alpha ion-ion streaming instability.  Include full kinetic electrons
        
        --D(omegahat) neglects the intrinsic Poisson response (i.e. assumes quasineutrality)
        --Assume adiabatic electron response
        --Also scales the dielectric function epsilon by a factor to yield D.
        --Here: assume omega is imaginary: and give the imaginary part (a real quantity)
            -->This last assumption means that D is real, so return a real number.

    Dimensionless variables:
    omegahat_i = imag(omega) / omega_pd
    kparallelhat = k_parallel * lambda_De
    kperphat = k_perp * lambda_De
    nu_d = n_d / n_e
    nu_alpha = n_alpha / n_e
    tau_e = T_e / T_d
    tau_alpha = T_alpha / T_d
    ui = V / w_Td
    """
    mu_e = 1 / 3672   # m_e / m_d
    mu_e2 = sqrt(2 / (nu_d * mu_e))
    tau_e2 = sqrt(2 / (nu_d * tau_e))
    tau_alpha2 = sqrt(tau_alpha / (nu_d * tau_e))
    khat = sqrt(kparallelhat^2 + kperphat^2)

    zeta_e = im*omegahat_i / (khat * mu_e2)
    zeta_d = im*omegahat_i / (khat * tau_e2)
    zeta_t = im*omegahat_i / (khat * sqrt(2/3) * tau_e2)
    zeta_aplus = (im*omegahat_i + kparallelhat * ui * tau_e2) / (khat * tau_alpha2)

    Zval_e = Z(zeta_e)
    Zval_d = Z(zeta_d)
    Zval_t = Z(zeta_t)
    Zval_aplus = Z(zeta_aplus)

    # scaled susceptibilities chi --- multiplied by k^2 lambda_Dd^2
    chi_e = 1 / (nu_d * tau_e) * (1 + real(zeta_e * Zval_e))
    chi_d = 1 + real(zeta_d * Zval_d)
    chi_t = 1 + real(zeta_t * Zval_t)
    chi_alpha = 4 * nu_alpha / (nu_d * tau_alpha) * 2 * real(1 + zeta_aplus * Zval_aplus)
    Dr = chi_e + chi_d + chi_t + chi_alpha
    return Dr
end


function D_aiisi(omegahat, kparallelhat, kperphat, nu_d, nu_alpha, tau_e, tau_alpha, ui)
    """dispersion function for alpha ion-ion streaming instability.
        --Include full kinetic electrons and do not assume quasineutrality
        --Also scales the dielectric function epsilon by a factor to yield D.

    Dimensionless variables:
    omegahat = omega / omega_pd
    kparallelhat = k_parallel * lambda_De
    kperphat = k_perp * lambda_De
    nu_d = n_d / n_e
    nu_alpha = n_alpha / n_e
    tau_e = T_e / T_d
    tau_alpha = T_alpha / T_d
    ui = V / w_Td
    """
    mu_e = 1 / 3672   # m_e / m_d
    mu_e2 = sqrt(2 / (nu_d * mu_e))
    tau_e2 = sqrt(2 / (nu_d * tau_e))
    tau_alpha2 = sqrt(tau_alpha / (nu_d * tau_e))
    khat = sqrt(kparallelhat^2 + kperphat^2)

    zeta_e = omegahat / (khat * mu_e2)
    zeta_d = omegahat / (khat * tau_e2)
    zeta_t = omegahat / (khat * sqrt(2/3) * tau_e2)
    zeta_aplus = (omegahat + kparallelhat * ui * tau_e2) / (khat * tau_alpha2)
    zeta_aminus = (omegahat - kparallelhat * ui * tau_e2) / (khat * tau_alpha2)

    Zval_e = Z(zeta_e)
    Zval_d = Z(zeta_d)
    Zval_t = Z(zeta_t)
    Zval_aplus = Z(zeta_aplus)
    Zval_aminus = Z(zeta_aminus)

    # scaled susceptibilities chi --- multiplied by k^2 lambda_Dd^2
    
    chi_e = 1 / (nu_d * tau_e) * (1 + zeta_e * Zval_e)
    chi_d = 1 + zeta_d * Zval_d
    chi_t = 1 + zeta_t * Zval_t
    chi_alpha = 4 * nu_alpha / (nu_d * tau_alpha) * (1 + zeta_aplus * Zval_aplus + 1 + zeta_aminus * Zval_aminus)
    D = khat^2 / (nu_d * tau_e)  +  chi_e + chi_d + chi_t + chi_alpha
    return D
end

function DDprime_aiisi(omegahat, kparallelhat, kperphat, nu_d, nu_alpha, tau_e, tau_alpha, ui)
    """dispersion function for alpha ion-ion streaming instability.
        --Include full kinetic electrons and do not assume quasineutrality
        --Also scales the dielectric function epsilon by a factor to yield D.
        --Return both D(omegahat), D'(omegahat)
    Dimensionless variables:
    omegahat = omega / omega_pd
    kparallelhat = k_parallel * lambda_De
    kperphat = k_perp * lambda_De
    nu_d = n_d / n_e
    nu_alpha = n_alpha / n_e
    tau_e = T_e / T_d
    tau_alpha = T_alpha / T_d
    ui = V / w_Td
    """
    mu_e = 1 / 3672   # m_e / m_d
    mu_e2 = sqrt(2 / (nu_d * mu_e))
    tau_e2 = sqrt(2 / (nu_d * tau_e))
    tau_alpha2 = sqrt(tau_alpha / (nu_d * tau_e))
    khat = sqrt(kparallelhat^2 + kperphat^2)

    zeta_e = omegahat / (khat * mu_e2)
    zeta_d = omegahat / (khat * tau_e2)
    zeta_t = omegahat / (khat * sqrt(2/3) * tau_e2)
    zeta_aplus = (omegahat + kparallelhat * ui * tau_e2) / (khat * tau_alpha2)
    zeta_aminus = (omegahat - kparallelhat * ui * tau_e2) / (khat * tau_alpha2)

    Zval_e, Zprimeval_e, Zhatval_e = Zprimehat(zeta_e)
    Zval_d, Zprimeval_d, Zhatval_d = Zprimehat(zeta_d)
    Zval_t, Zprimeval_t, Zhatval_t = Zprimehat(zeta_t)
    Zval_aplus, Zprimeval_aplus, Zhatval_aplus = Zprimehat(zeta_aplus)
    Zval_aminus, Zprimeval_aminus, Zhatval_aminus = Zprimehat(zeta_aminus)
    

    # scaled susceptibilities chi --- multiplied by k^2 lambda_Dd^2
    chi_e = 1 / (nu_d * tau_e) * (1 + zeta_e * Zval_e)
    chiprime_e = 1 / (nu_d * tau_e * khat * mu_e2) * Zhatval_e
    
    chi_d = 1 + zeta_d * Zval_d
    chiprime_d = 1 / (khat * tau_e2)  * Zhatval_d
    
    chi_t = 1 + zeta_t * Zval_t
    chiprime_t = 1 / (khat * sqrt(2/3) * tau_e2)  * Zhatval_t
    
    chi_alpha = 4 * nu_alpha / (nu_d * tau_alpha) * (1 + zeta_aplus * Zval_aplus + 1 + zeta_aminus * Zval_aminus)
    chiprime_alpha = 4 * nu_alpha / (nu_d * tau_alpha * khat * tau_alpha2) * (Zhatval_aplus + Zhatval_aminus)

    D = khat^2 / (nu_d * tau_e)  +  chi_e + chi_d + chi_t + chi_alpha
    Dprime = chiprime_e + chiprime_d + chiprime_t + chiprime_alpha
    return D, Dprime
end


# function D_ionacoustic(omega_n, kn, tau)
#     """Dispersion function for ion acoustic waves.
    
#     evaluate D(omega_n) for adiabatic electrons, Maxwellian ions.
    
#     omega_n = omega/omega_pi, k_n = k*lambda_de,  tau = Ti/Z*Te 

#     Also evaluate D'(omega_n) = dD/domega_n
    
#     Add inputs: dimensionless parameters"""
#     dxi_i = 1 / (kn * sqrt(an.ArbReal(2) * tau))   # calculate d(xi_i)/d(omega_n)
#     xi_i = omega_n * dxi_i
#     Zval, Zprimeval, Zhatval = Zprimehat(xi_i)
#     D = (kn^2 + 1) * tau  +  1  +  xi_i * Zval
#     Dprime = Zhatval * dxi_i
#     return D, Dprime
# end

# function D_foote(omega_n, eta, ui, tauinv)
#     """Use dispersion function in Foote & Kulsrud for ion-ion streaming instability.
        
#         --D(omega_n) in Foote neglects the intrinsic Poisson response (i.e. assumes quasineutrality)
#         --Also scales the dielectric function epsilon by a factor to yield D.
    
#     Here, evaluate both D(omega_n) and D'(omega_n) for adiabatic electrons, Maxwellian ions.
    
#     omega_n = omega/omega_pi, eta = k*V/omega_pi,  ui=w_ti/V,  tauinv = Ti/Te 
    
#     Add inputs: dimensionless parameters"""
#     dxi_i = 1 / (eta*ui)  # calculate d(xi_i)/d(omega_n)
#     xi_iplus = (omega_n + eta)  * dxi_i
#     xi_iminus = (omega_n - eta)  * dxi_i
    
#     Zval_plus, Zprimeval_plus, Zhatval_plus = Zprimehat(xi_iplus)
#     Zval_minus, Zprimeval_minus, Zhatval_minus = Zprimehat(xi_iminus)
#     D = tauinv +  1  +  (xi_iplus * Zval_plus  +  xi_iminus * Zval_minus) / 2
#     Dprime = dxi_i * (Zhatval_plus + Zhatval_minus) / 2 
#     return D, Dprime
# end

# function D_issi(omegahat, kparallelhat, kperphat, ui, tau, mu)
#     """
#     Use dispersion relation for ion-ion streaming instability (ISSI).  Do not make the 2 assumptions in
#     Foote and Kulrsud of quasineutrality or Boltzmann electron response.

#     Here, evaluate both D(omegahat) and D'(omegahat)

#     Dimensionless variables:
#     omegahat = omega / omega_pi
#     kparallelhat = k_parallel * lambda_De
#     kperphat = k_perp * lambda_De
#     ui = U / w_Ti
#     tau = T_e / T_i
#     mu = m_i / m_e
#     """
#     mu2 = 2 * sqrt(mu)
#     tau2 = 2 / sqrt(tau)
#     khat = sqrt(kparallelhat^2 + kperphat^2)
#     dzeta_i = 1 / (khat * tau2)  # calculate d(xi_i)/d(omega_hat)
#     dzeta_e = 1 / (khat * mu2)  # calculate d(xi_e)/d(omega_hat)

#     zeta_e = omegahat * dzeta_e
#     zeta_iplus = (omegahat + kparallelhat * ui * tau2) * dzeta_i
#     zeta_iminus = (omegahat - kparallelhat * ui * tau2) * dzeta_i

#     Zval_e, Zprimeval_e, Zhatval_e = Zprimehat(zeta_e)
#     Zval_iplus, Zprimeval_iplus, Zhatval_iplus = Zprimehat(zeta_iplus)
#     Zval_iminus, Zprimeval_iminus, Zhatval_iminus = Zprimehat(zeta_iminus)

#     D = khat^2 - Zprimeval_e/2  - tau/4 * (Zprimeval_iminus + Zprimeval_iplus)
#     Dprime = Zhatval_e * dzeta_e  + tau/2 * dzeta_i * (Zhatval_iminus + Zhatval_iplus)

#     return D, Dprime
# end

# function D_ion_acoustic_instability(omegahat, khat, tau, ui)
#     Z = an.ArbReal(1)
#     mime = an.ArbReal(1836)
#     dxi_i = sqrt(tau/2) / khat  # calculate d(xi_i)/d(omega_hat)
#     dxi_e = sqrt(Z/mime) / khat # calculate d(xi_e)/d(omega_hat)

#     xi_e = (omegahat - khat * ui) * dxi_e
#     xi_i = omegahat * dxi_i

#     Zval_e, Zprimeval_e, Zhatval_e = Zprimehat(xi_e)
#     Zval_i, Zprimeval_i, Zhatval_i = Zprimehat(xi_i)

#     D = khat^2 / tau  +  (1 + xi_e * Zval_e) / tau  + 1  +  xi_i * Zval_i
#     Dprime = dxi_e * Zhatval_e / tau   +   dxi_i * Zhatval_i
#     return D, Dprime
# end

end