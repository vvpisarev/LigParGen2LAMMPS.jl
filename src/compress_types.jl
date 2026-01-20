function compress_types!(mol::LigParGenMolecule, typenames::Vector{Symbol})
    type_param = Tuple{Symbol, Float64, Float64, Float64}[]

    for (name, m, (eps, sig, t)) in zip(mol.masses, typenames, mol.pair_coeffs)
        k = findfirst(==((name, m, eps, sig),), type_param)
        if k === nothing
            push!(type_param, (name, m, eps, sig))
            k = length(type_param)
        end

        for (itype, type) in pairs(mol.types)
            if type == t
                mol.types[itype] = k
            end
        end
    end

    resize!(mol.masses, length(type_param))
    resize!(mol.pair_coeffs, length(type_param))

    for (t, (m, eps, sig)) in pairs(type_param)
        mol.masses[t] = m
        mol.pair_coeffs[t] = (eps, sig, t)
    end

    return mol
end

function remove_null_dihedrals!(mol::LigParGenMolecule)
    null_dtypes = Int[]
    for k in keys(mol.dihed_coeffs)
        coeffs = mol.dihed_coeffs[k]
        if all(iszero, coeffs)
            push!(null_dtypes, k)
        end
    end

    deleteat!(mol.diheds, null_dtypes)
    deleteat!(mol.dihed_coeffs, null_dtypes)

    null_itypes = Int[]
    for k in keys(mol.improper_coeffs)
        coeffs = mol.improper_coeffs[k]
        if coeffs[1] == 0
            push!(null_itypes, k)
        end
    end

    deleteat!(mol.improps, null_itypes)
    deleteat!(mol.improper_coeffs, null_itypes)
    return mol
end
