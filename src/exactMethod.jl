

function ExactLin(expr_fct::Ef,x1::Real,x2::Real, e::ErrorType; bounding = Best() ::BoundingType, 
                    ConcavityChanges = [Inf]::Array{Float64,1} )
    if x1 >= x2
        return Float64[]
    end
   
    pwl = Array{LinearPiece}(undef, 0)

    if false
        f1 = expr_fct(x1)
        f2 = expr_fct(x2)
        if max(abs(x1), abs(f1)) < 1e-7
            xt = x2
            ft = expr_fct(xt)
            while abs(ft) > 1e-2 * abs(f2)
                xt = (x1 + xt) / 2
                ft = expr_fct(xt)
            end
            lp = LinearPiece(x1, xt, 0, 0, x -> 0.0)
            println("first piece approximated to zero in interval [$x1, $xt]")
            push!(pwl, lp)
            x1 = xt
            # return vcat(lp, ExactLin(expr_fct, xt, x2, e; bounding = bounding, ConcavityChanges = ConcavityChanges))
        end
    end

    if ConcavityChanges == [Inf]
        ConcavityChanges = ConcavitySplit(x1,x2,expr_fct)
    end
    ConcavityChanges = [x1;ConcavityChanges;x2] 
    ConcavityChanges = sort(unique(ConcavityChanges)) # make sure that the bounds are there
    cor = CorridorFromInfo(x1, x2,expr_fct, e, bounding)[3:end-1]
    
    x2Temp = -1
    i=1
    while x1 < x2
        
        #find next concavity change
        i = searchsortedfirst(ConcavityChanges,x1)
        ConcavityChanges[i] == x1 ? x2Temp = ConcavityChanges[i+1] : x2Temp = ConcavityChanges[i]
        
        #first apply the algortihm for the convex/concave segement
        pwl = [pwl;Linearize(expr_fct, x1, x2Temp, e,bounding = bounding,ConcavityChanges=[] )] 
        pwl[end].xMax == x2 && break

        #Replace the last segment by one obtain with the general algorithm
        pwl[end] = ExactPiece(pwl[end].xMin,x2,cor...)

        pwl[end].xMax == x2 && break

        x1 = pwl[end].xMax
    end
    return pwl
end