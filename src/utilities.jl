
"""
isContinuous(pwl, ε = 1e-5)

Determine whether a pwl function is continuous up to a numerical precision of ε.
# Arguments
- `plw` : pwl function
- `ε` : numerical precision used to detect if the end endpoints of two segments are equal (to detect discontinuities)

"""
function isContinuous(pwl, ε = EPS) 

    
    for i in 1:length(pwl)-1
        
        temp = pwl[i].xMax
        
        if pwl[i](temp) - pwl[i+1](temp) > ε
            return false
        end
        
    end
    
    true
    
end





"""
breakpoints(pwl, ε = 1e-5)

Outputs the breakpoints needed for several solvers (CPLEX, Gurobi,...) to natively model PWL functions. 
# Arguments
- `plw` : pwl function
- `ε` : numerical precision used to detect if the end endpoints of two segments are equal (to detect discontinuities)

"""
function breakpoints(pwl, ε = EPS) 

    
    bpx = [pwl[1].xMin]

    
    for i in 1:length(pwl)-1
        
        temp = pwl[i].xMax
        push!(bpx,temp)
        
        if pwl[i](temp) - pwl[i+1](temp) > ε
            push!(bpx , temp)
        end
        
    end
    
    push!(bpx,pwl[end].xMax)
    
    return bpx
    
    
end

function find_zeros(f::Ef, x1::Real, x2::Real)
    rts = roots(f, interval(x1, x2))
    if isempty(rts) return Float64[]
    else return [z.region.bareinterval.lo for z in rts]
    end
end

function find_zero(f::Ef, x1::Real, x2::Real)
    return find_zeros(f, x1, x2)[begin]
end

function maximize(f::Ef, x1::Real, x2::Real)::ScalarOptResult
    df(x) = Derive(f)(x)
    zs = find_zeros(x -> df(x), x1, x2)
    append!(zs, [x1, x2])
    return ScalarOptResult(argmax(x -> f(x), zs), maximum(x -> f(x), zs))
end

function minimize(f::Ef, x1::Real, x2::Real)
    res = maximize(x -> -f(x), x1, x2)
    return ScalarOptResult(res.x, -res.val)
end

function get_scale(f::Ef, x1::Real, x2::Real)::Float64
    resmax = maximize(f, x1, x2)
    resmin = minimize(f, x1, x2)
    vmax = resmax.val
    vmin = resmin.val
    return max(abs(vmax), abs(vmin))
end

function scale_function(f::Ef, s::Real, x1::Real, x2::Real)::Ef
    if f isa Expr
        return :( eval(f) * s )
    elseif f isa Function
        return y -> f(x1 + y * (x2 - x1)) / s
    end
end

function invert_function(f::Ef)::Ef
    if f isa Expr
        return :( - eval(f) )
    elseif f isa Function
        return x -> - f(x)
    end 
end

function is_mostly_negative(f::Ef, x1::Real, x2::Real)::Bool
    res = maximize(f, x1, x2)
    return res.val < EPS
end

function unscale_linearpiece(lp::LinearPiece, s::Real, x1::Real, x2::Real)::LinearPiece
    ymin = x1 + lp.xMin * (x2 - x1) 
    ymax = x1 + lp.xMax * (x2 - x1)
    ap = s * lp.a / (x2 - x1)
    bp = s * (lp.b - lp.a * x1 / (x2 - x1))
    return LinearPiece(ymin, ymax, ap, bp, x -> ap * x + bp)
end

function construct_constant_piece(f::Ef, x1::Real, x2::Real, bounding::BoundingType)
    if bounding == Under()
        b = minimize(f, x1, x2).val
    elseif bounding == Over()
        b = maximize(f, x1, x2).val
    else
        b = (minimize(f, x1, x2).val + maximize(f, x1, x2).val) / 2
    end
    return LinearPiece(x1, x2, 0.0, b, x -> b)
end
    
function invert_linearpiece(lp::LinearPiece, inv::Int64)
    if inv == 1 return lp
    else
        return LinearPiece(lp.xMin, lp.xMax, -lp.a, -lp.b, x -> -lp.fct(x))
    end
end

function reduce_infeasibilities(f::Ef, lp::LinearPiece, bounding::BoundingType)::LinearPiece
    if bounding isa Best
        return lp
    elseif bounding isa Under
        res = minimize(x -> f(x) - lp(x), lp.xMin, lp.xMax)
        val = res.val
        if val < -EPS
            println("translating piece by $(val)")
            newlp = lp + val
            return newlp
        else
            return lp
        end
    else
        res = minimize(x -> lp(x) - f(x), lp.xMin, lp.xMax)
        val = res.val
        if val < -EPS
            println("translating piece by $(-val)")
            newlp = lp - val
            return newlp
        else
            return lp
        end
    end
end

function (pwl::Vector{LinA.LinearPiece})(x::Real, bounding, eps = EPS)::LinearizationEval
    if x < pwl[1].xMin - eps || x > pwl[end].xMax + eps
        throw(DomainError(x, "argument must be in the domain of the function"))
    end
    f, l = 1, length(pwl)
    m = 0
    while f <= l
        m = floor(Int64, (f + l) / 2)
        p = pwl[m]
        if x >= p.xMin - eps && x <= p.xMax + eps
            if bounding isa Best
                return LinearizationEval(p(x), m)
            else
                break
            end
        elseif x < p.xMin - 1e-9
            l = m - 1
        else
            f = m + 1
        end
    end
    optval = bounding isa Under ? 1e+20 : -1e+20
    optfn = bounding isa Under ? min : max
    isless_fn(u, v)  = bounding isa Under ? u < v : u > v
    optm = m
    for i = m : -1 : 1
        p = pwl[i]
        if x < p.xMin - eps || x > p.xMax + eps
            break
        elseif isless(p(x), optval)
            optval = p(x)
            optm = i
        end
    end
    for i = m + 1 : length(pwl)
        p = pwl[i]
        if x < p.xMin - eps || x > p.xMax + eps
            break
        elseif isless(p(x), optval)
            optval = p(x)
            optm = i
        end
    end
    return LinearizationEval(optval, optm)
end

function compute_limits_at_zero(g::Ef, x1::Real, x2::Real)::Vector{IntervalArithmetic.Interval}
    eps = 1e-4
    zs = IntervalRootFinding.roots(g, interval(x1, x2))
    rts = [(z.region.bareinterval.lo + z.region.bareinterval.hi) / 2 for z in zs]
    zints = IntervalArithmetic.Interval[]
    l(x) = g(x) - eps
    u(x) = g(x) + eps
    for (i, r) in enumerate(rts)
        x0 = i > 1 ? rts[i - 1] + eps : x1
        xf = i < length(rts) ? rts[i + 1] - eps : x2
        zl = IntervalRootFinding.roots(x -> l(x), interval(x0, xf))
        zr = IntervalRootFinding.roots(x -> u(x), interval(x0, xf))
        ul = isempty(zl) ? typemax(Float64) : zl[begin].region.bareinterval.lo
        ur = isempty(zr) ? typemax(Float64) : zr[begin].region.bareinterval.lo
        @assert(min(ul, ur) < typemax(Float64), "cannot be both roots empty! $zl, $zr")
        if max(ul, ur) < typemax(Float64)
            y0 = min(ul, ur) 
            yf = max(ul, ur)
        else
            u = min(ul, ur)
            if u > r
                y0 = r
                yf = u
            else
                y0 = u
                yf = r
            end
        end
        y0 = max(x1, min(y0, r - eps))
        yf = min(x2, max(yf, r + eps))
        push!(zints, interval(y0, yf))
    end
    sort!(zints; lt = (u, v) -> isstrictless(u, v))
    return zints
end
function compute_nonzero_intervals(g::Ef, x1::Real, x2::Real, zints::Vector{IntervalArithmetic.Interval})::Vector{IntervalArithmetic.Interval}
    left = interval(x1, x2)
    nzints = IntervalArithmetic.Interval[]
    for zint in zints
        chop = interiordiff(left, zint)
        if !isatomic(chop[begin])
            push!(nzints, chop[begin])
        end
        left = chop[end]
    end
    if !isatomic(left)
        push!(nzints, left)
    end
    return nzints
end