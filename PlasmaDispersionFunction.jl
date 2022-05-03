"""
Solving plasma dispersion relations and determining stability.  For:
    --Modified two-stream instability (MSTI)

Calculating zeros of a dispersion function D(omega) or dielectric function epsilon(omega)
* Solve for local zeros using an iterative root finder (e.g., fixed-point or Newton-like)
* For analytic D, find zeros within a certain contour with the Davies method

Other methods of determining stability
* Count zeros of analytic epsilon within a contour using the argument principle (or use real line + small imag part)
    --This is essentially the Nyquist diagram
"""

module PlasmaDispersionFunction
import SpecialFunctions



function Z(z)
    """Plasma dispersion function Z(z).
    Faddeeva function is w(z) = erfcx(-iz).
    Plasma dispersion function is Z(z) = i * sqrt(pi) * w(z)

    This is valid at large |z|, Im(z) > 0
    """
    out = im * sqrt(pi) * SpecialFunctions.erfcx(-im*z)
    return out
end

function Zprime(z)
    """More robust version of Z'(z).  Get better precision at large |z| by using the asymptotic series.
    For imag(z) >= 0, use the series at |z| >= 1e2.
    """
    y = imag(z)
    if y >= 0 && abs(z) >= 100
        if y > 0
            out1 = 0
        else
            sigma = 1
            out1 = im * sqrt(pi) * sigma * -2 * z * exp(-z^2)
        end
        out2 = 1/z^2 + 3/2 / z^4 + 15/4 / z^6 + 105/8 / z^8
        return out1 + out2
    else
        return -2 * (1 + z*Z(z))
    end
end

function Z2(z)
    """Robust version.  Return Z(z), Z'(z), and Z''(z) all together.  (Note Z'' = -2*Zhat)
    Take care at large |z| for Im(z) > 0"""
    Zval = Z(z)  # this evaluation is good for large |z|, Im(z) > 0

    # for Z', Zhat, use asymptotic series at Im(z) > 0 and |z| > 100
    y = imag(z)
    if y >= 0 && abs(z) > 100
        if y > 0
            Zprime_out1 = 0
            Zhat_out1 = 0
        else
            sigma = 1
            Zprime_out1 = im * sqrt(pi) * sigma * -2 * z * exp(-z^2)
            out1zZprime = im * sqrt(pi) * sigma * -2 * z^2 * exp(-z^2)
            Z_out1 = im * sqrt(pi) * sigma * exp(-z^2)
            Zhat_out1 = z*Zprime_out1 + Z_out1
        end
        Zprime_out2 = 1/z^2 + 3/2 / z^4 + 15/4 / z^6 + 105/8 / z^8
        Zprimeval = Zprime_out1 + Zprime_out2
        Zhat_out2 = 1/z^3 * (1 + 3/z^2 + 45/4 / z^4)
        Zhatval = Zhat_out1 + Zhat_out2
    else
        Zprimeval = -2 * (1 + z * Zval)
        Zhatval = z * Zprimeval + Zval
    end
    Zprime2 = -2 * Zhatval  # Zprime2 = Z''.
    return Zval, Zprimeval, Zprime2
end

function Zprimehat(z)
    """Robust version.  Return Z(z), Zprime(z), and Zhat(z) all together.
    Take care at large |z| for Im(z) > 0"""
    Zval = Z(z)  # this evaluation is good for large |z|, Im(z) > 0

    # for Z', Zhat, use asymptotic series at Im(z) > 0 and |z| > 100
    y = imag(z)
    if y >= 0 && abs(z) > 100
        if y > 0
            Zprime_out1 = 0
            Zhat_out1 = 0
        else
            sigma = 1
            Zprime_out1 = im * sqrt(pi) * sigma * -2 * z * exp(-z^2)
            out1zZprime = im * sqrt(pi) * sigma * -2 * z^2 * exp(-z^2)
            Z_out1 = im * sqrt(pi) * sigma * exp(-z^2)
            Zhat_out1 = z*Zprime_out1 + Z_out1
        end
        Zprime_out2 = 1/z^2 + 3/2 / z^4 + 15/4 / z^6 + 105/8 / z^8
        Zprimeval = Zprime_out1 + Zprime_out2
        Zhat_out2 = 1/z^3 * (1 + 3/z^2 + 45/4 / z^4)
        Zhatval = Zhat_out1 + Zhat_out2
    else
        Zprimeval = -2 * (1 + z * Zval)
        Zhatval = z * Zprimeval + Zval
    end
    return Zval, Zprimeval, Zhatval
end

# -----------------------------------------------------------
#        Reference functions (likely not used in practice)
# ----------------------------------------------------------


function Zasymptotic(z)
    """Asymptotic series evaluation of Z(z).  Not needed since Z works well for Im(z) >= 0 even for large |z|.
        Here as reference."""
    y = imag(z)
    if y > 0
        sigma = 0
        out1 = 0
    elseif y == 0
        sigma = 1
        out1 = im * sqrt(pi) * sigma * exp(-z^2)
    else
        sigma = 2
        out1 = im * sqrt(pi) * sigma * exp(-z^2)
    end
    out2 = -1/z * (1 + 1/2 / z^2 + 3/4 / z^4 + 15/8 / z^6 + 105/16 / z^8)
    return out1 + out2
end

function Zprime_basic(z)
    """Derivatve Z'(z).  Note that at large |z| there can be precision loss, so multiple precision helps"""
    out = -2 * (1 + z * Z(z))
    return out
end



function Zhat_basic(z)
    """Derivative Zhat(z) = d(zZ)/dz = zZ' + Z.  

    **Note that Zhat(z) = -Z''(z) / 2
    At large |z|, there are 2 rounds of precision loss (one in Zprimeval, one in Zhat), so large loss of precision."""
    Zval = Z(z)
    Zprimeval = -2 * (1 + z * Zval)
    out = z * Zprimeval + Zval
    return out
end

function Zhat(z)
    """More robust version of Z'(z).  Get better precision at large |z| by using the asymptotic series.
    For imag(z) >= 0, use the series at |z| >= 1e2.
    """
    y = imag(z)
    if y >= 0 && abs(z) >= 100
        if y > 0
            out1 = 0
        else
            sigma = 1
            out1zZprime = im * sqrt(pi) * sigma * -2 * z^2 * exp(-z^2)
            out1Z = im * sqrt(pi) * sigma * exp(-z^2)
            out1 = out1zZprime + out1Z
        end
        out2 = 1/z^3 * (1 + 3/z^2 + 45/4 / z^4)
    else
        return z*Zprime(z) + Z(z)
    end
end

function Zprimehat_basic(z)
    """Return Z(z), Zprime(z), and Zhat(z) all together.  Not being careful at large |z|"""
    Zval = Z(z)
    Zprimeval = -2 * (1 + z * Zval)
    Zhatval = z * Zprimeval + Zval
    return Zval, Zprimeval, Zhatval
end

end