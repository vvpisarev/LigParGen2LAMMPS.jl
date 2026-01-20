# Original contribution by Ivan M. Varaksin
# Refactoring by Vasily V. Pisarev

"""
    switch_ff!(mol::LigParGenMolecule, typenames, ff_type::Symbol)

Switch LJ, bond, angle and dihedral parameters to the desired force field.

Known force fields: :opls_aa, :opls_aa_2020
"""
function switch_ff!(mol::LigParGenMolecule, typenames, ff_type::Symbol)
    if !(ff_type in (:opls_aa, :opls_aa_2020, :opls_aa_ligpargen))
        error("Unrecognized FF type, allowed types: :opls_aa, :opls_aa_2020")
    end

    if ff_type == :opls_aa
        switch_ff_to_oplsaa!(mol, typenames)
    elseif ff_type == :opls_aa_2020
        switch_ff_to_opls2020!(mol, typenames)
    end
    return mol
end

function switch_ff_to_oplsaa!(mol::LigParGenMolecule, typenames)
    ct_types = (:CT0, :CT1, :CT2, :CT3, :CT4)
    for k in eachindex(mol.pair_coeffs)
        typeind = mol.pair_coeffs[k][3]
        if typenames[typeind] == :HC
            mol.pair_coeffs[k] = (0.030, 2.5, typeind)
        elseif typenames[typeind] in ct_types
            mol.pair_coeffs[k] = (0.066, 3.5, typeind)
        end
    end
    for k in eachindex(mol.dihed_coeffs)
        d = mol.diheds[k]
        if all(ind -> typenames[ind] in ct_types, d[1])
            mol.dihed_coeffs[k] = (1.3, -0.05, 0.2, 0.0)
        end
    end
    return mol
end

function switch_ff_to_opls2020!(mol::LigParGenMolecule, typenames)
    ct_types = (:CT0, :CT1, :CT2, :CT3, :CT4)
    for k in keys(mol.pair_coeffs)
        typeind = mol.pair_coeffs[k][3]
        if typenames[typeind] == :HC
            mol.pair_coeffs[k] = (0.026, 2.48, typeind)
        elseif typenames[typeind] == :CT4
            mol.pair_coeffs[k] = (0.060, 3.57, typeind)
        elseif typenames[typeind] == :CT3
            mol.pair_coeffs[k] = (0.066, 3.55, typeind)
        elseif typenames[typeind] == :CT2
            mol.pair_coeffs[k] = (0.066, 3.51, typeind)
        elseif typenames[typeind] == :CT1
            mol.pair_coeffs[k] = (0.066, 3.5, typeind)
        elseif typenames[typeind] == :CT0
            mol.pair_coeffs[k] = (0.066, 3.5, typeind)
        end
    end
    for k in keys(mol.dihed_coeffs)
        d = mol.diheds[k]
        dihed_types = getindex.(Ref(typenames), d)
        if dihed_types == (:CT3, :CT2, :CT2, :CT3) # n-butane
            mol.dihed_coeffs[k] = (1.10, -0.20, 0.2, 0.0)
        elseif all(in(ct_types), dihed_types)
            mol.dihed_coeffs[k] = (0.85, -0.20, 0.2, 0.0)
        end
    end
    return mol
end

function switch_charge!(mol::LigParGenMolecule, typenames, charge_model::Symbol)
    charge_model == :cm1a && return mol
    charge_model in (:opls_aa, :cm1a) || error("Unknown charge model: $(charge_model)")
    for (k, name) in pairs(typenames)
        if name == :HC
            mol.charges[k] = 0.06
        elseif name == :CT4
            mol.charges[k] = -0.24
        elseif name == :CT3
            mol.charges[k] = -0.18
        elseif name == :CT2
            mol.charges[k] = -0.12
        elseif name == :CT1
            mol.charges[k] = -0.06
        elseif name == :CT0
            mol.charges[k] = 0.0
        end
    end
    return mol
end
