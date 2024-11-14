
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

function scale_linearpiece(lp::LinearPiece, s::Real, x1::Real, x2::Real)::LinearPiece
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
