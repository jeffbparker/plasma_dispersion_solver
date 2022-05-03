"""
Tools for doing numerical continuation when finding zeros of functions as a function of parameter.

Continuation in one parameter is possible.
"""

module continuation

import Roots
import Polynomials

function polyfiteval(x_in, y_in, x_out)
    """Fit a polynomial of degree N-1 to the N points in (x_in, y_in).  Then, evaluate on x_out.
    
    Inputs:
        x_in     input points, x values (N-element vector)
        y_in     input points, y values (N-element vector)
        x_out    output point, x value (scalar)
    Outputs:
        y_out    output point, y value (scalar)
    """
    # fit N-1 degree polynomial to the N points
    f = Polynomials.fit(x_in, y_in)
    y_out = f(x_out)
    return y_out
end

function get_prior_points(vec, j, Nmax)
    """Get up to Nmax points preceding index j from input vector.
    
    We wish to select up to Nmax prior points (immediately preceding index j).
    If j <= Nmax, fewer than Nmax points are available, but select as many points as possible.

    Inputs:
        vec        input vector (1d array)
        j          current index (integer)
    Outputs:
        vec_prior  prior points (1d array, up to size Nmax)
    """
    if j <= Nmax  # fewer than Nmax points available
        vec_prior = vec[1:j-1]
    elseif j >= Nmax + 1 # at least Nmax prior points available
        vec_prior = vec[j-Nmax:j-1]
    end
    return vec_prior
end

function do_continuation_vectorinput(bvec, xsol, f, fprime)
    """continuation using input vector bvec rather than start value, stop value, and number of points.
    xsol is the solution at bvec[1].
    
    Here, f and fprime are functions with signature fun(x, b) where x is the complex variable (or vector)
    and b is the parameter.  If the original problem has multiple parameters, f and fprime must be parameterized
    to match this signature.
    
    Note that xsolvec[1] is not solved for---it is just xsol.

    Inputs:
        bvec        array of parameter values b for which to solve the continuation problem.  (1d array)
        xsol        Value x at which f(xsol, b[1]) = 0.  This must be a good value; it is not checked it is a solution here (scalar, real or complex)
        f           function for which to find zeros, signature f(x, b)
        fprime      f' = partial f/partial x, signature fprime(x, b) 
    
    Outputs:
        xsolvec     array of solution values x, for which f(xsolvec[j], bvec[j]) = 0  (1d array)
    """
    N = length(bvec)
    xsolvec = zeros(typeof(xsol), N)
    xsolvec[1] = xsol
    for j in range(2, N)
        #get 3 previous points from bvec, xvec
        Nmax = 3
        inputs_prior = get_prior_points(bvec, j, Nmax)
        outputs_prior = get_prior_points(xsolvec, j, Nmax)

        # use previous points to fit and find guess for next solution value
        inputval_next = bvec[j]
        outputval_next_guess = polyfiteval(inputs_prior, outputs_prior, inputval_next)

        # solve fzero with guess
        g(x) = f(x, inputval_next)
        gprime(x) = fprime(x, inputval_next)
        output_soln = try
            Roots.find_zero((g,gprime), outputval_next_guess, Roots.Newton())
        catch e
            println("Couldn't converge somewhere; ending early")
            return xsolvec
        end
        xsolvec[j] = output_soln
    end
    return xsolvec
end

function do_continuation(bstart, bstop, Nb, xstart, f, fprime)
    """continuation using input vector bvec rather than start value, stop value, and number of points.
    xsol is the solution at bvec[1].
    
    Here, f and fprime are functions with signature fun(x, b) where x is the complex variable (or vector)
    and b is the parameter.  If the original problem has multiple parameters, f and fprime must be parameterized
    to match this signature.
    
    Note that xsolvec[1] is not solved for---it is just xsol.

    Inputs:
        bstart      start value b at which to solve the continuation problem.  (scalar))
        bstop       stop value of b for continuation (scalar)
        Nb          Number of points b, inclusive of bstart, bstop (integer)
        xstart      Value x at which f(xstart, bstart) = 0.  This must be a valid value; it is not checked it is a solution here (scalar, real or complex)
        f           function for which to find zeros, signature f(x, b)
        fprime      f' = partial f/partial x, signature fprime(x, b) 
    
    Outputs:
        bvec        array of parameter values b, at which continuation problem is solved.  (1d array)
        xsolvec     array of solution values x, for which f(xsolvec[j], bvec[j]) = 0  (1d array)
    """
    bvec = LinRange(bstart, bstop, Nb)
    xsolvec = do_continuation_vectorinput(bvec, xstart, f, fprime)
    # xsolvec = zeros(typeof(xstart), Nb)
    # xsolvec[1] = xstart
    # for j in range(2, N)
    #     #get 3 previous points from bvec, xvec
    #     Nmax = 3
    #     inputs_prior = get_prior_points(bvec, j, Nmax)
    #     outputs_prior = get_prior_points(xsolvec, j, Nmax)

    #     # use previous points to fit and find guess for next solution value
    #     inputval_next = bvec[j]
    #     outputval_next_guess = polyfiteval(inputs_prior, outputs_prior, inputval_next)

    #     # solve fzero with guess
    #     g(x) = f(x, inputval_next)
    #     gprime(x) = fprime(x, inputval_next)
    #     output_soln = try
    #         Roots.find_zero((g,gprime), outputval_next_guess, Roots.Newton())
    #     catch e
    #         println("Couldn't converge somewhere; ending early")
    #         return xsolvec
    #     end
    #     xsolvec[j] = output_soln
    # end
    return bvec, xsolvec
end

function do_continuation_leftright(bstart, xstart, bleft, Nleft, bright, Nright, f, fprime)
    """Continuation to both left and right from some middle point (xstart, bstart) which should
    satisfy f(xstart, bstart)=0.

    Inputs:
        bstart      start value b at which to solve the continuation problem.  (scalar)
        xstart      Value x at which f(xstart, bstart) = 0 (this is NOT checked). (scalar, real or complex)
        bleft       stop value of b for continuation to the left (scalar)
        Nleft       Number of points for continuation to the left, from bstart to bleft (integer)
        bright      stop value of b for continuation to the right (scalar)
        Nright      Number of points for continuation to the right, from bstart to bright (integer)
        f           function for which to find zeros, signature f(x, b)
        fprime      f' = partial f/partial x, signature fprime(x, b)
    """
    # Perform continuation to the left and to the right
    bvec_left, xsoln_left = do_continuation(bstart, bleft, Nleft, xstart, f, fprime)
    bvec_right, xsoln_right = do_continuation(bstart, bright, Nright, xstart, f, fprime)
    
    # reverse bvec_left, xsoln_left so they are in increasing order
    bvec_left = reverse(bvec_left)
    xsoln_left = reverse(xsoln_left)

    # concatenate the two solutions, deleting the redundant point
    bvec = vcat(bvec_left[1:end-1], bvec_right)
    xsoln = vcat(xsoln_left[1:end-1], xsoln_right)
    return bvec, xsoln
end

end