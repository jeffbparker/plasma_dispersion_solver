"""
Tools for using the argument principle to count and find zeros.
    Basic argument principle
    Compute the zeros using Delves & Lynesse; Davies
    Squircle contour from QuaLiKiz

This module provides a path for doing these computations in double arithmetic.
"""
module tools
import Polynomials

function centered_diff(u, dx)
    """Compute du/dx.
      du/dx is computed using centered differences on the same grid as u.  For the edge points, one-point differences are used.
    
    Inputs:
      u         profile (array)
      dx        grid spacing (scalar)
    
    Outputs:
      dudx      (array, same length as u)
    """
    dudx = zero(u)
    dudx[1] = (u[2] - u[1]) / dx
    dudx[2:end-1] = (u[3:end] - u[1:end-2]) / (2*dx)
    dudx[end] = (u[end] - u[end-1]) / dx
    return dudx
end

function numerical_complex_derivative(f, a, h)
    """Compute df/dz at a, using 4-point Cauchy integral formula with circle of radius h.
    Good default value: h = sqrt(eps)
    
    Inputs:
        f           function: C to C (Complex numbers)
        a           value at which to evaluate f'
        h           radius of circle to use in Cauchy integral formula
    Outputs:
        fprime      approximation to f'(a) (complex number)
    """
    z1 = a + h
    z2 = a - h
    z3 = a + h*im
    z4 = a - h*im
    fprime = (f(z1) - f(z2) -im * (f(z3) - f(z4))) / (4*h) 
    return fprime
end

function linspace(start, stop, N)
    """N points between start and stop, inclusive.
    Should work with float or ArbNumerics.  
    
    Presumably the built-in function LinRange does the same thing.  For some reason
    I was unsatisfied using it for ArbNumerics.ArbReal, but LinRange seems to work just fine.
    """
    dx = (stop - start) / (N - 1)
    x = zeros(typeof(dx), N)
    x[1] = start
    for j in 2:N
        x[j] = x[j-1] + dx
    end
    return x
end

function theta_circle(N)
    """Return N points theta, uniformly distributed on the circle in [0, 2pi)."""
    arb2pi = 2 * pi
    theta = zeros(Float64, N)
    for j in 1:N
        theta[j] = arb2pi * (j-1) / N
    end
    return theta
end

# later: figure out how to do an abstract type, to use Float64 or ArbReal
struct Sqparams
    """Squircle parameters for specifying the shifted squircle contour."""
    a::Float64
    Rx::Float64
    Ry::Float64
    rx::Float64
    epsilon::Float64
    smallnum::Float64
end

function sq(theta, params::Sqparams)
    """Squircle map for input on the unit circle.
    For scalar input theta, not array.

    Inputs:
        theta      angle on unit circle (scalar)
        params
    Outputs:
        xa         x-coordinate of base (unshifted) squircle, (scalar)
        ya         y-coordinate of base (unshifted) squircle, (scalar)
    """
    a, smallnum = params.a, params.smallnum
    c = cos(theta)
    s = sin(theta)
        
    # handle values specially where sin(theta) or cos(theta) is very small to avoid divide by zero.
    arbpi = pi
    theta2 = mod(theta, arbpi / 2)
    if abs(theta2) < smallnum || abs(theta2 - arbpi/2) < smallnum
        xa = a * c
        ya = a * s
        return xa, ya
    end
    prefactor = sign(c*s) / sqrt(2)
    h = sqrt(1 - sqrt(1 - 4 * a*a * c*c * s*s))
    factor = prefactor * h

    xa = factor / s
    ya = factor / c
    return xa, ya
end

function dsq_dtheta(theta, params::Sqparams)
    """Derivative of squircle map for input on the unit circle.
    If the squircle map is (xa, ya) = gamma(theta), this returns (dxa/dtheta, dya/dtheta)
    For scalar input, not vector.
    
    Inputs:
        theta       angle on unit circle (scalar)
        params
    Outputs:
        dxa_dtheta  x-coordinate of derivative of base (unshifted) squircle, (scalar)
        dya_dtheta  y-coordinate of derivative of base (unshifted) squircle, (scalar)
    """
    a, smallnum = params.a, params.smallnum
    
    c = cos(theta)
    s = sin(theta)
    
    # handle values specially where sin(theta) or cos(theta) is very small to avoid divide by zero
    arbpi = pi
    theta2 = mod(theta, arbpi / 2)
    if abs(theta2) < smallnum || abs(theta2 - arbpi/2) < smallnum
        xa = -a * s
        ya = a * c
        return xa, ya
    end

    prefactor = sign(c*s) / sqrt(2)
    temp = sqrt(1 - 4 * a*a * c*c * s*s)
    h = sqrt(1 - temp)
    dh_dtheta = 2*a*a* s * c * (c*c - s*s) / (h*temp)
    factor = prefactor * h
    dxa_dtheta = prefactor * (dh_dtheta * s - h * c) / (s*s)
    dya_dtheta = prefactor * (dh_dtheta * c + h * s) / (c*c)
    
    return dxa_dtheta, dya_dtheta
end

function xyprime(theta, params::Sqparams)
    """Compute x', y' of the shifted squircle

    Inputs:
        theta      angle on unit circle (scalar)
        params
    Outputs:
        xp         x-coordinate of shifted squircle, (scalar)
        yp         y-coordinate of shifted squircle, (scalar)
    """
    a, Rx, Ry, rx, epsilon = params.a, params.Rx, params.Ry, params.rx, params.epsilon
    xa, ya = sq(theta, params)
    xp = Rx + (rx/a) * xa
    yp = Ry + (Ry - epsilon)/a * ya
    return xp, yp
end

function dxyprime_dtheta(theta, params::Sqparams)
    """Compute dx'/dtheta, dy'/dtheta of the shifted squircle
    
    Inputs:
        theta           angle on unit circle (scalar)
        params
    Outputs:
        dxp_dtheta      x-coordinate of shifted squircle, (scalar)
        dyp_dtheta      y-coordinate of shifted squircle, (scalar)
    """
    a, Rx, Ry, rx, epsilon = params.a, params.Rx, params.Ry, params.rx, params.epsilon
    dxa_dtheta, dya_dtheta = dsq_dtheta(theta, params)
    dxp_dtheta = (rx/a) * dxa_dtheta
    dyp_dtheta = (Ry - epsilon)/a * dya_dtheta
    return dxp_dtheta, dyp_dtheta
end

function gamma(theta, params::Sqparams)
    """Shifted squircle.  Output is complex number z=x+iy rather than tuple (x,y)
    
    Inputs:
        theta      angle on unit circle (scalar)
        params
    Outputs:
        z          coordinate of shifted squircle in complex form x + iy, (complex scalar)
    """
    xp, yp = xyprime(theta, params)
    z = xp + im*yp
    return z
end

function dgamma_dtheta(theta, params::Sqparams)
    """Shifted squircle, derivative.  Output is complex number dz/dtheta rather than (dx/dtheta, dy/dtheta)
    
    Inputs:
        theta       angle on unit circle (scalar)
        params
    Outputs:
        dz_dtheta   derivative of shifted squircle in complex form, (complex scalar)
    """
    dxp_dtheta, dyp_dtheta = dxyprime_dtheta(theta, params)
    dz_dtheta = dxp_dtheta + im*dyp_dtheta
    return dz_dtheta
end

# ----------------------- Functions for counting zeros --------------------------- #
function count_zeros_from_vals(vals)
    """Using the argument principle, count zeros of an analytic function inside a contour using values
    evaluated on the contour.  Assume no poles inside the contour
    
    Inputs:
        vals            values of function evaluated on contour (array, complex numbers)
    Outputs:
        winding_number  number of zeros detected inside (integer)
    """
    args = angle.(vals)   # argument of complex number computed with atan2
    d_args = diff(args)
    ind_plus = (d_args .> 1.5*pi)
    ind_minus = (d_args .< -1.5*pi)
    winding_number_plus = -count(ind_plus)   # count(ind) <--> np.count_nonzero(ind)
    winding_number_minus = count(ind_minus)

    winding_number = winding_number_plus + winding_number_minus
    
    # check that all of the rest of the jumps are less than pi/2
    ind2 = (abs.(d_args) .<= 1.5*pi)
    d_args_small = d_args[ind2]
    
    ind3 = abs.(d_args_small) .<= 0.5*pi
    tempcount = count(ind3)
    
    @assert tempcount == length(d_args_small) "Jump in argument too large detected.  Try increasing number of sampling points."
    return winding_number
end

function count_zeros_inside_sq_contour(fun, N, params::Sqparams)
    """Using the argument principle, count zeros of an analytic function fun inside a contour using values
    evaluated on the contour.  Assume no poles inside the contour
    
    Inputs:
        fun             complex analytic function within squircle contour f(z), takes complex number
        N               Number of points to evalute function on squircle contour (integer)
        params          squircle params
    Outputs:
        num_zeros       number of zeros detected inside the squircle contour (integer)
    """
    theta = theta_circle(N)
    vals = zeros(ComplexF64, N)
    for j = 1:N
        point_on_squircle = gamma(theta[j], params)
        vals[j] = fun(point_on_squircle)
    end

    num_zeros = count_zeros_from_vals(vals)
    return num_zeros
end


# ------------------ Functions for integrating on the squircle contour ----------- #
function integrate_fun_on_squircle(fun, N, params::Sqparams)
    """Integrate a generic function on the squircle.  Here for reference.
    
    Integral is carried out using uniformly spaced points on unit circle, mapped to squircle,
    with integration performed via the trapezoidal method (converging exponentially fast).
    
    Inputs:
        fun             complex analytic function within squircle contour f(z), takes complex number
        N               Number of points to use on squircle contour (integer)
        params          squircle params
    Outputs:
        result          approximation to value of integral of fun on the squircle contour (complex scalar)
    """
    theta = theta_circle(N)
    dtheta = (2 * pi) / N
    
    fout = zeros(ComplexF64, N)
    for j in 1:N
        fout[j] = fun(gamma(theta[j], params)) * dgamma_dtheta(theta[j], params)
    end
    result = dtheta * sum(fout)
    return result
end

function compute_Sarray(f, fprime, num_zeros, Npts, params::Sqparams)
    """Compute all of the Sn simultaneously for Davies method.  The number of zeros of f must have been found already.
    Integrate multiple functions at the same time so as to reuse function evaluations for efficiency.
        
    Inputs:
        f             complex analytic function f(z) for which to find zeros, takes complex number
        fprime        complex function f'(z) = df/dz, takes complex number
        num_zeros     number of zeros of f within the squircle contour (integer)
        Npts          Number of points to use on squircle contour for integration  (integer)
        params        squircle params
    Outputs:
        Sn            approximation to value of integrals Sn on the squircle contour (array of complex numbrs)    
    """
    theta = theta_circle(Npts)
    Sn = zeros(ComplexF64, num_zeros)    
    for j in 1:Npts
        z = gamma(theta[j], params)
        dg_dtheta = dgamma_dtheta(theta[j], params)
        fval = f(z)
        fprimeval = fprime(z)
        s1_integrand = z * fprimeval / fval * dg_dtheta   # integrand for S1 is z fprime(z)/f 
        Sn[1] += s1_integrand
        s_integrand = s1_integrand
        for n in 2:num_zeros
            sn_integrand = z * s_integrand  # integrand for Sn is z^n fprime(z)/f
            Sn[n] += sn_integrand
        end
    end
    Sn = Sn ./ (Npts * im)
    return Sn
end

function compute_Sarray(f_fprime, num_zeros, Npts, params::Sqparams)
    """Compute all of the Sn simultaneously for Davies method.  The number of zeros of f must have been found already.
    Integrate multiple functions at the same time so as to reuse function evaluations for efficiency.
      
    Here, f_fprime is a function of z that computes and returns both f(z), f'(z) [computing these together may be efficient]

    Inputs:
        f_fprime      function of z that computes and returns both f(z), f'(z).  Type of z should be complex number
        num_zeros     number of zeros of f within the squircle contour (integer)
        Npts          Number of points to use on squircle contour for integration  (integer)
        params        squircle params
    Outputs:
        Sn            approximation to value of integrals Sn on the squircle contour (array of complex numbers)    
    """
    theta = theta_circle(Npts)
    Sn = zeros(ComplexF64, num_zeros)    
    for j in 1:Npts
        z = gamma(theta[j], params)
        dg_dtheta = dgamma_dtheta(theta[j], params)
        fval, fprimeval = f_fprime(z)
        s1_integrand = z * fprimeval / fval * dg_dtheta # integrand for S1 is z^n fprime(z)/f
        Sn[1] += s1_integrand
        s_integrand = s1_integrand
        for n in 2:num_zeros
            sn_integrand = z * s_integrand  # integrand for Sn is z^n fprime(z)/f
            Sn[n] += sn_integrand
        end
    end
    Sn = Sn ./ (Npts * im)
    return Sn
end

# ------------ Functions for finding zeros using the Davies polynomial method --------- #
function Sarray_to_Aarray(Sarray)
    """Transform Sn coefficients to the An coefficients using the recursive method.
    (Newton's identities for sums of powers in polynomial roots)
    
    Inputs:
        Sarray =     [S1, S2, ..., SN]     (complex array)
    
    Outputs:
        Aarray = [A0, A1, A2, ..., AN]     (complex array)
    """
    num_zeros = length(Sarray)
    Aarray = zeros(ComplexF64, num_zeros + 1)
    Aarray[1] = 1   # A0 = 1
    for k in 1:num_zeros
        # Compute A_k
        tempsum = 0
        for j in 1:k
            tempsum += Aarray[k-j + 1] * Sarray[j]
        end
        Aarray[k+1] = -tempsum / k
    end
    return Aarray
end

function roots_of_poly(coeffarray)
    """Find roots of polynomial using the Polynomial package.

    Inputs:
        coeffarray   polynomial coefficients in usual order---lowest degree first. (complex array)
    Outputs:
        roots        roots of polynomial (complex array)
    """
    p = Polynomials.Polynomial(coeffarray)
    roots = Polynomials.roots(p)
    return roots
end

function root_with_largest_imag_part(roots)
    """Given roots object returned by the Polynomial package, return the root with largest imaginary part
    
    Inputs:
        roots        roots of polynomial. (complex array)
    Outputs:
        root         root of polynomial with largest imaginary part (scalar Complex)
    """
    # sort by imaginary part (descending order)
    roots = sort(roots, by=imag, rev=true)
    root = roots[1]
    return root
end

function max_growth_rate(f, fprime, num_zeros, Npts, params::Sqparams)
    """Use Davies method on a dispersion relation f(omega).
    Get the root with maximum growth rate (largest imaginary party omega_i)
    Returns as ArbReals (omega_r, omega_i).

    Inputs:
        f             complex analytic function f(z) for which to find zeros, takes complex number
        fprime        complex function f'(z) = df/dz, takes complex number
        num_zeros     number of zeros of f within the squircle contour (integer)
        Npts          Number of points to use on squircle contour for integration  (integer)
        params        squircle params
    Outputs:
        omega_r       Real part of eigenmode (real frequency)  (scalar)
        omega_i       Imaginary part of eigenmode (growth rate)  (scalar)
    """

    Sarray = compute_Sarray(f, fprime, num_zeros, Npts, params)

    Aarray = Sarray_to_Aarray(Sarray)
    coeffArray = reverse(Aarray)
    roots = roots_of_poly(coeffArray)
    root = root_with_largest_imag_part(roots)   # returns a BigFloat
    omega_r = real(root)
    omega_i = imag(root)
    return omega_r, omega_i
end

function max_growth_rate(f_fprime, num_zeros, Npts, params::Sqparams)
    """Use Davies method on a dispersion relation f(omega).
    Get the root with maximum growth rate (largest imaginary party omega_i)
    Returns as ArbReals (omega_r, omega_i).  Note that BigFloat is used in intermediate, for finding polynomial roots.

    Here, f_fprime is a function of z that computes and returns both f(z), f'(z)
    Inputs:
        f_fprime      function of z that computes and returns both f(z), f'(z).  Type of z should be complex number
        num_zeros     number of zeros of f within the squircle contour (integer)
        Npts          Number of points to use on squircle contour for integration  (integer)
        params        squircle params
    Outputs:
        omega_r       Real part of eigenmode (real frequency)  (scalar)
        omega_i       Imaginary part of eigenmode (growth rate)  (scalar)
    """
    Sarray = compute_Sarray(f_fprime, num_zeros, Npts, params)

    Aarray = Sarray_to_Aarray(Sarray)
    coeffArray = reverse(Aarray)
    roots = roots_of_poly(coeffArray)
    root = root_with_largest_imag_part(roots)
    omega_r = real(root)
    omega_i = imag(root)
    return omega_r, omega_i
end

function all_roots(f_fprime, num_zeros, Npts, params::Sqparams)
    """Use Davies method on a dispersion relation f(omega).
    Get the root with maximum growth rate (largest imaginary party omega_i)
    Returns as ArbReals (omega_r, omega_i).  Note that BigFloat is used in intermediate, for finding polynomial roots.

    Here, f_fprime is a function of z that computes and returns both f(z), f'(z)
    Inputs:
        f_fprime      function of z that computes and returns both f(z), f'(z).  Type of z should be complex number
        num_zeros     number of zeros of f within the squircle contour (integer)
        Npts          Number of points to use on squircle contour for integration  (integer)
        params        squircle params
    Outputs:
        omega_r       Real part of eigenmode (real frequency)  (scalar)
        omega_i       Imaginary part of eigenmode (growth rate)  (scalar)
    """
    Sarray = compute_Sarray(f_fprime, num_zeros, Npts, params)
    Aarray = Sarray_to_Aarray(Sarray)
    coeffArray = reverse(Aarray)
    roots = roots_of_poly(coeffArray)
    return roots
end

end