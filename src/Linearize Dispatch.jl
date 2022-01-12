include("Heuristic.jl")
include("exact method.jl")


"""
    Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType, ConcavityChanges = [Inf]::Array{Float64,1})

Makes an optimal piecewise Linear approximation of expr_fct from x1 to x2. The result will be an array of `LinearPiece`. Note that the array is directly callable as a function.
# Arguments
- `expr_fct` : function to linearize (either an expresion or a native julia function)
- `x1` : from
- `x2` : to
- `e` : error type either Absolute() or Relative()

# Optional Arguments
- `bounding` : `Under()` for an underestimation, `Over()` for an overestimation, `Best()` for estimation that can go under or over the function, `UnderOver()` for both an underestimation and an overestimation in two separate estimations ( saves a lot of overhead over calling `Over()` and `Under()` separately). By default, it uses `Best()`
- `ConcavityChanges` : Concavity changes in the function. If not given, they will be computed automatically which, in rare cases, can lead to precision errors if the concavity is asymptotic to zero. 

# Example
```julia-repl
julia> pwl = Linearize(:(x^2),0,2,Absolute(0.1))
3-element Vector{LinA.LinearPiece}:
 0.894427190999916 x -0.1 from 0.0 to 0.894427190999916
 2.683281572999748 x -1.7000000000000006 from 0.894427190999916 to 1.7888543819998326
 4.736067977499794 x -5.372135954999589 from 1.7888543819998326 to 2.0

julia> pwl(1)
0.9832815729997475
```
"""
function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType,
ConcavityChanges = [Inf]::Array{Float64,1})

    return HeuristicLin(expr_fct,x1,x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges)
end

function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm::HeuristicLin; bounding = Best() ::BoundingType,
    ConcavityChanges = [Inf]::Array{Float64,1})

        return HeuristicLin(expr_fct,x1,x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges)
end

function Linearize(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType,algorithm::ExactLin; bounding = Best() ::BoundingType,
    ConcavityChanges = [Inf]::Array{Float64,1})

        return exactLin(expr_fct,x1,x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges)
end
