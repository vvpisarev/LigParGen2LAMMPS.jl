const COMBINATIONS = let
    sums = (1, -1, 2, -2, 3, -3, 4, -4, 5, -5)
    ab_pairs = [(a, b) for s in sums for a in 0:sign(s):s for b in (s-a, a-s)]
    unique!(ab_pairs)
    sort!(ab_pairs; by=((a, b),)->max(abs(a), abs(b)))
end

function balance_charges!(
    mol::LigParGenMolecule, desired_charge, dq=1e-4; charge_diff_thresh
)
    charge = mol.charges
    pcoeff = mol.pair_coeffs
    type = mol.types
    charge_groups = Vector{Int}[]

    charge_vals = Float64[]

    for (i, q) in pairs(charge)
        k = findfirst(abs(qval - q) < charge_diff_thresh for qval in charge_vals)
        if k === nothing
            push!(charge_vals, q)
            push!(charge_groups, [i])
        else
            inds = charge_groups[k]
            type_i = type[i]
            eps, sig = pcoeff[type_i]
            epsg, sigg = pcoeff[type[first(inds)]]
            if (eps, sig) == (epsg, sigg)
                push!(inds, i)
                gr = @view charge[inds]
                charge_vals[k] = sum(gr) / length(gr)
            else
                push!(charge_vals, q)
                push!(charge_groups, [i])
            end
        end
    end

    for group in charge_groups
        grview = @view charge[group]
        ave_charge = sum(grview) / length(grview)
        grview .= dq * round(ave_charge / dq)
    end

    charge_correction = float(desired_charge) - sum(charge)

    nq = round(Int, charge_correction / dq)

    nq == 0 && return mol

    @info "" nq string(charge_groups)

    for (a, b) in COMBINATIONS
        for i1 in eachindex(charge_groups)
            gr1 = charge_groups[i1]
            len1 = length(gr1)
            for i2 in nextind(charge_groups, i1):lastindex(charge_groups)
                gr2 = charge_groups[i2]
                len2 = length(gr2)
                if a * len1 + b * len2 == nq
                    @views charge[gr1] .+= a * dq
                    @views charge[gr2] .+= b * dq
                    @info "" sum(charge)
                    return mol
                end
            end
        end
    end

    best_group = argmax(length(gr) for gr in charge_groups)
    gr = charge_groups[best_group]

    if length(gr) > abs(nq)
        for i in range(div(length(gr) - abs(nq), 2) + 1; length=abs(nq))
            @views charge[gr][i] += dq * sign(nq)
        end
    else
        @views charge[gr] .+= nq / length(gr) * dq
    end

    @info "" sum(charge)
    return mol
end
