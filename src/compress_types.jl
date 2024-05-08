function compress_types!(mol::Molecule)
    type_param = NTuple{3,Float64}[]

    for (m, (eps, sig, t)) in zip(mol.masses, mol.pair_coeffs)
        k = findfirst(==((m, eps, sig),), type_param)
        if k === nothing
            push!(type_param, (m, eps, sig))
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
