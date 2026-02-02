"""
    export_mol(f, mol::Molecule[; ff, charge_model])

Write the molecular topology of `mol` into `f`. `f` can be an I/O stream or a filename.
"""
function export_mol(io::IO, mol::Molecule; ff=:opls_aa_ligpargen, charge_model=:cm1a)
    alt_base = deepcopy(mol.base)
    remove_null_dihedrals!(alt_base)
    switch_ff!(alt_base, mol.typenames, ff)
    switch_charge!(alt_base, mol.typenames, charge_model)
    __write_mol(io, alt_base)
    # mol.compress_types && compress_types!(alt_base, mol.typenames)
    # mol.compress_btypes && compress_btypes!(alt_base)
    # mol.compress_atypes && compress_atypes!(alt_base)
    # mol.compress_dtypes && compress_dtypes!(alt_base)
    # mol.compress_itypes && compress_itypes!(alt_base)
    return nothing
end

function export_mol(
    fname::AbstractString, mol::Molecule,
    ;
    ff=:opls_aa_ligpargen, charge_model=:cm1a,
)
    open(fname, "w") do io
        export_mol(io, mol)
    end
end

function __write_mol(io::IO, mol::LigParGenMolecule)
    println(io, mol.comment)
    println(io, length(mol.coords), " atoms")
    println(io, length(mol.bonds), " bonds")
    println(io, length(mol.angles), " angles")
    println(io, length(mol.diheds), " dihedrals")
    println(io, length(mol.improps), " impropers")

    print_coords(io, mol.coords)
    print_types(io, mol.types)
    print_masses(io, mol.masses, mol.types)
    print_charges(io, mol.charges)

    !isempty(mol.bonds) && print_bonds(io, mol.bonds)
    !isempty(mol.angles) && print_angles(io, mol.angles)
    !isempty(mol.diheds) && print_dihedrals(io, mol.diheds)
    !isempty(mol.improps) && print_impropers(io, mol.improps)
end

"""
    export_mol(f, ligpargenfile::AbstractString)

Read the molecular topology from file named `ligpargenfile` and write it into `f`. `f` can
    be an I/O stream or a filename.
"""
function export_mol(f, mol::AbstractString; ff=:opls_aa_ligpargen, charge_model=:cm1a, kw...)
    export_mol(f, read_lpg_data(mol, kw...); ff, charge_model)
end

"""
    export_ff(f, mol::Molecule; include_charge=true, ff=nothing, charge=nothing)

Write the forcefield parameters of `mol` into `f`. `f` can be an I/O stream or a filename.
    `ff` can be specified as `:opls_aa` or `:opls_aa_2020`, `charge` can be specified as
    `:opls_aa`.
"""
function export_ff(f::IO, mol::Molecule; ff=:opls_aa_ligpargen, charge_model=:cm1a)
    alt_base = deepcopy(mol.base)
    remove_null_dihedrals!(alt_base)
    switch_ff!(alt_base, mol.typenames, ff)
    switch_charge!(alt_base, mol.typenames, charge_model)
    __write_ff(f, alt_base)
    return nothing
end

function __write_ff(io, mol::LigParGenMolecule)
    println(io,
    """
    pair_style        lj/cut/coul/long 12.0
    bond_style        harmonic
    angle_style       harmonic
    """ *
    (isempty(mol.diheds) ? "" : "dihedral_style    opls\n") *
    (isempty(mol.improps) ? "" : "improper_style    cvff\n") *
    """
    special_bonds     lj/coul 0.0 0.0 0.5
    pair_modify       mix geometric
    """)

    for (id, m) in pairs(mol.masses)
        join(io, ("mass          ", id, m, '\n'), ' ')
    end

    !isempty(mol.pair_coeffs) && println(io)

    for type in keys(mol.pair_coeffs)
        ε, σ = mol.pair_coeffs[type]
        join(io, ("pair_coeff    ", type, type, ε, σ, '\n'), ' ')
    end

    !isempty(mol.bond_coeffs) && println(io)

    for btype in keys(mol.bond_coeffs)
        k, r0 = mol.bond_coeffs[btype]
        join(io, ("bond_coeff    ", btype, k, r0, '\n'), ' ')
    end

    !isempty(mol.angle_coeffs) && println(io)

    for atype in 1:length(mol.angle_coeffs)
        k, θ0 = mol.angle_coeffs[atype]
        join(io, ("angle_coeff    ", atype, k, θ0, '\n'), ' ')
    end

    !isempty(mol.dihed_coeffs) && println(io)

    for dtype in 1:length(mol.dihed_coeffs)
        c1, c2, c3, c4 = mol.dihed_coeffs[dtype]
        join(io, ("dihedral_coeff    ", dtype, c1, c2, c3, c4, '\n'), ' ')
    end

    !isempty(mol.improper_coeffs) && println(io)

    for itype in 1:length(mol.improper_coeffs)
        k, d, n = mol.improper_coeffs[itype]
        join(io, ("improper_coeff    ", itype, k, d, n, '\n'), ' ')
    end

    println(io)

    for (k, q) in pairs(mol.charges)
        join(io, ("set    type ", k, " charge ", round(q; digits=5), '\n'))
    end
end

function export_ff(fname::AbstractString, mol::Molecule; kw...)
    open(fname, "w") do io
        export_ff(io, mol; kw...)
    end
end

"""
    export_ff(f, ligpargenfile::AbstractString)

Read the molecular topology from file named `ligpargenfile` and write its forcefield
    into `f`. `f` can be an I/O stream or a filename.
"""
function export_ff(f, mol::AbstractString; ff=:opls_aa_ligpargen, charge_model=:cm1a, kw...)
    export_ff(f, read_lpg_data(mol; kw...); ff, charge_model)
end

"""
    export_mol_and_ff(fmask::AbstractString, mol::Molecule)

Write files "fmask.txt" and "fmask.ff" with molecular topology
and forcefield parameters of `mol`, respectively.
"""
function export_mol_and_ff(
    fmask::AbstractString, mol::Molecule,
    ;
    ff=:opls_aa_ligpargen, charge_model=:cm1a,
)
    alt_base = deepcopy(mol.base)
    remove_null_dihedrals!(alt_base)
    switch_ff!(alt_base, mol.typenames, ff)
    switch_charge!(alt_base, mol.typenames, charge_model)
    @sync begin
        @async __write_ff("$fmask.ff", mol)
        @async __write_mol("$fmask.txt", mol)
    end
    return
end

"""
    export_mol_and_ff(fmask::AbstractString, ligpargenfile::AbstractString)

Read the molecular topology and forcefield parameters from `ligpargenfile` and write them
    into files "fmask.txt" and "fmask.ff", respectively.
"""
export_mol_and_ff(fmask, mol::AbstractString; kw...) = export_mol_and_ff(fmask, read_lpg_data(mol; kw...))

function print_coords(io, coords::Vector)
    println(io, "\nCoords\n")
    for i in 1:length(coords)
        x, y, z = coords[i]
        join(io, (i, x, y, z, '\n'), ' ')
    end
end

function print_types(io, types::Vector)
    println(io, "\nTypes\n")
    for i in 1:length(types)
        join(io, (i, types[i], '\n'), ' ')
    end
end

function print_charges(io, charges::Vector)
    println(io, "\nCharges\n")
    for i in 1:length(charges)
        join(io, (i, round(charges[i]; digits=14), '\n'), ' ')
    end
end

function print_masses(io, masses::Vector, types::Vector)
    println(io, "\nMasses\n")
    for (i, t) in pairs(types)
        join(io, (i, masses[t], '\n'), ' ')
    end
end

function print_bonds(io, bonds::Vector)
    println(io, "\nBonds\n")
    for i in keys(bonds)
        (a1, a2) = bonds[i]
        type = i
        join(io, (i, type, a1, a2, '\n'), ' ')
    end
end

function print_angles(io, angles::Vector)
    println(io, "\nAngles\n")
    for i in keys(angles)
        (a1, a2, a3) = angles[i]
        type = i
        join(io, (i, type, a1, a2, a3, '\n'), ' ')
    end
end

function print_dihedrals(io, dihed::Vector)
    println(io, "\nDihedrals\n")
    for i in keys(dihed)
        (a1, a2, a3, a4) = dihed[i]
        type = i
        join(io, (i, type, a1, a2, a3, a4, '\n'), ' ')
    end
end

function print_impropers(io, impropers::Vector)
    println(io, "\nImpropers\n")
    for i in keys(impropers)
        (a1, a2, a3, a4) = impropers[i]
        type = i
        join(io, (i, type, a1, a2, a3, a4, '\n'), ' ')
    end
end

"""
    export_data(f, mol::Molecule[; ff, charge_model])

Write the molecular topology and forcefield parameters of `mol` as LAMMPS data file.
"""
function export_data(io::IO, mol::Molecule; ff=:opls_aa_ligpargen, charge_model=:cm1a)
    alt_base = deepcopy(mol.base)
    remove_null_dihedrals!(alt_base)
    switch_ff!(alt_base, mol.typenames, ff)
    switch_charge!(alt_base, mol.typenames, charge_model)
    __write_data(io, alt_base)
    return nothing
end

function export_data(path::AbstractString, mol::Molecule; kw...)
    open(path, "w") do io
        export_data(io, mol; kw...)
    end
end

function __write_data(io::IO, mol::LigParGenMolecule)
    println(io, mol.comment)
    println(io)
    println(io, length(mol.coords), " atoms")
    println(io, length(mol.bonds), " bonds")
    println(io, length(mol.angles), " angles")
    isempty(mol.diheds) || println(io, length(mol.diheds), " dihedrals")
    isempty(mol.improps) || println(io, length(mol.improps), " impropers")

    println(io)

    println(io, length(mol.pair_coeffs), " atom types")
    println(io, length(mol.bond_coeffs), " bond types")
    println(io, length(mol.angle_coeffs), " angle types")
    isempty(mol.dihed_coeffs) || println(io, length(mol.diheds), " dihedral types")
    isempty(mol.improper_coeffs) || println(io, length(mol.improps), " improper types")

    println(io)

    println(io, mol.box_info)

    println(io, "Masses\n")

    for (i, t) in pairs(mol.types)
        join(io, (i, mol.masses[t], '\n'), ' ')
    end

    println(io, "\nPair Coeffs # lj/cut/coul/long\n")
    for type in keys(mol.pair_coeffs)
        ε, σ = mol.pair_coeffs[type]
        join(io, (type, ε, σ, '\n'), ' ')
    end

    if !isempty(mol.bond_coeffs)
        println(io, "\nBond Coeffs # harmonic\n")

        for btype in keys(mol.bond_coeffs)
            k, r0 = mol.bond_coeffs[btype]
            join(io, (btype, k, r0, '\n'), ' ')
        end
    end

    if !isempty(mol.angle_coeffs)
        println(io, "\nAngle Coeffs # harmonic\n")

        for atype in 1:length(mol.angle_coeffs)
            k, θ0 = mol.angle_coeffs[atype]
            join(io, (atype, k, θ0, '\n'), ' ')
        end
    end

    if !isempty(mol.dihed_coeffs)
        println(io, "\nDihedral Coeffs # opls\n")

        for dtype in 1:length(mol.dihed_coeffs)
            c1, c2, c3, c4 = mol.dihed_coeffs[dtype]
            join(io, (dtype, c1, c2, c3, c4, '\n'), ' ')
        end
    end

    if !isempty(mol.improper_coeffs)
        println(io, "\nImproper Coeffs # cvff\n")
        for itype in 1:length(mol.improper_coeffs)
            k, d, n = mol.improper_coeffs[itype]
            join(io, (itype, k, d, n, '\n'), ' ')
        end
    end

    println(io, "\nAtoms\n")
    for (i, (x, y, z)) in pairs(mol.coords)
        join(io, (i, 1, mol.types[i], round(mol.charges[i]; digits=6), x, y, z, '\n'), ' ')
    end

    println(io, "\nBonds\n")
    for (i, (a, b)) in pairs(mol.bonds)
        join(io, (i, i, a, b, '\n'), ' ')
    end

    println(io, "\nAngles\n")
    for (i, (a, b, c)) in pairs(mol.angles)
        join(io, (i, i, a, b, c, '\n'), ' ')
    end

    if !isempty(mol.diheds)
        println(io, "\nDihedrals\n")
        for (i, (a, b, c, d)) in pairs(mol.diheds)
            join(io, (i, i, a, b, c, d, '\n'), ' ')
        end
    end

    if !isempty(mol.improps)
        println(io, "\nImpropers\n")
        for (i, (a, b, c, d)) in pairs(mol.improps)
            join(io, (i, i, a, b, c, d, '\n'), ' ')
        end
    end
end
