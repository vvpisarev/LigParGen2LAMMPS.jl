export write_mol, write_ff, write_mol_and_ff

"""
    write_mol(f, mol::Molecule)

Write the molecular topology of `mol` into `f`. `f` can be an I/O stream or a filename.
"""
function write_mol(io, mol::Molecule)
    println(io, mol.comment)
    println(io, length(mol.coords), " atoms")
    println(io, length(mol.bonds), " bonds")
    println(io, length(mol.angles), " angles")
    println(io, length(mol.diheds), " dihedrals")
    println(io, length(mol.improps), " impropers")

    print_coords(io, mol.coords)
    print_types(io, mol.types)
    print_masses(io, mol.masses)
    print_charges(io, mol.charges)

    !isempty(mol.bonds) && print_bonds(io, mol.bonds)
    !isempty(mol.angles) && print_angles(io, mol.angles)
    !isempty(mol.diheds) && print_dihedrals(io, mol.diheds)
    !isempty(mol.improps) && print_impropers(io, mol.improps)
end

function write_mol(fname::AbstractString, mol::Molecule)
    open(fname, "w") do io
        write_mol(io, mol)
    end
end

"""
    write_mol(f, ligpargenfile::AbstractString)

Read the molecular topology from file named `ligpargenfile` and write it into `f`. `f` can
    be an I/O stream or a filename.
"""
write_mol(f, mol::AbstractString) = write_mol(f, read_lpg_data(mol))

"""
    write_ff(f, mol::Molecule)

Write the forcefield parameters of `mol` into `f`. `f` can be an I/O stream or a filename.
"""
function write_ff(io, mol::Molecule)
    println(io,
    """
    pair_style        lj/cut/coul/long 12.0
    bond_style        harmonic
    angle_style       harmonic
    dihedral_style    opls
    improper_style    cvff
    special_bonds     lj/coul 0.0 0.0 0.5
    pair_modify       mix geometric
    """)

    for type in 1:length(mol.pair_coeffs)
        itype = findlast(==(type), mol.types)
        m = mol.masses[itype]
        join(io, ("mass          ", type, m, '\n'), ' ')
    end

    !isempty(mol.pair_coeffs) && println(io)

    for type in 1:length(mol.pair_coeffs)
        ε, σ = mol.pair_coeffs[type]
        join(io, ("pair_coeff    ", type, type, ε, σ, '\n'), ' ')
    end

    !isempty(mol.bond_coeffs) && println(io)

    for btype in 1:length(mol.bond_coeffs)
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
end

function write_ff(fname::AbstractString, mol::Molecule)
    open(fname, "w") do io
        write_ff(io, mol)
    end
end

"""
    write_ff(f, ligpargenfile::AbstractString)

Read the molecular topology from file named `ligpargenfile` and write its forcefield
    into `f`. `f` can be an I/O stream or a filename.
"""
write_ff(f, mol::AbstractString) = write_ff(f, read_lpg_data(mol))

"""
    write_mol_and_ff(fmask::AbstractString, mol::Molecule)

Write files "fmask.txt" and "fmask.ff" with molecular topology
and forcefield parameters of `mol`, respectively.
"""
function write_mol_and_ff(fmask::AbstractString, mol::Molecule)
    @sync begin
        @async write_ff("$fmask.ff", mol)
        @async write_mol("$fmask.txt", mol)
    end
    return
end

"""
    write_mol_and_ff(fmask::AbstractString, ligpargenfile::AbstractString)

Read the molecular topology and forcefield parameters from `ligpargenfile` and write them
    into files "fmask.txt" and "fmask.ff", respectively.
"""
write_mol_and_ff(fmask, mol::AbstractString) = write_mol_and_ff(fmask, read_lpg_data(mol))

function print_coords(io, coords)
    println(io, "\nCoords\n")
    for i in 1:length(coords)
        x, y, z = coords[i]
        join(io, (i, x, y, z, '\n'), ' ')
    end
end

function print_types(io, types)
    println(io, "\nTypes\n")
    for i in 1:length(types)
        join(io, (i, types[i], '\n'), ' ')
    end
end

function print_charges(io, charges)
    println(io, "\nCharges\n")
    for i in 1:length(charges)
        join(io, (i, charges[i], '\n'), ' ')
    end
end

function print_masses(io, masses)
    println(io, "\nMasses\n")
    for i in 1:length(masses)
        join(io, (i, masses[i], '\n'), ' ')
    end
end

function print_bonds(io, bonds)
    println(io, "\nBonds\n")
    for i in 1:length(bonds)
        (a1, a2), type = bonds[i]
        join(io, (i, type, a1, a2, '\n'), ' ')
    end
end

function print_angles(io, angles)
    println(io, "\nAngles\n")
    for i in 1:length(angles)
        (a1, a2, a3), type = angles[i]
        join(io, (i, type, a1, a2, a3, '\n'), ' ')
    end
end

function print_dihedrals(io, dihed)
    println(io, "\nDihedrals\n")
    for i in 1:length(dihed)
        (a1, a2, a3, a4), type = dihed[i]
        join(io, (i, type, a1, a2, a3, a4, '\n'), ' ')
    end
end

function print_impropers(io, impropers)
    println(io, "\nImpropers\n")
    for i in 1:length(impropers)
        (a1, a2, a3, a4), type = impropers[i]
        join(io, (i, type, a1, a2, a3, a4, '\n'), ' ')
    end
end
