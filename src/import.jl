Base.@kwdef struct Molecule
    coords::Vector{NTuple{3,Float64}} = NTuple{3,Float64}[]
    types::Vector{Int16} = Int16[]
    typenames::Vector{Symbol} = Symbol[]
    charges::Vector{Float64} = Float64[]
    masses::Vector{Float64} = Float64[]
    bonds::Vector{Pair{NTuple{2,Int16},Int16}} = Pair{NTuple{2,Int16},Int16}[] #bond is (i, j) => btype
    angles::Vector{Pair{NTuple{3,Int16},Int16}} = Pair{NTuple{3,Int16},Int16}[] #angle is (i, j, k) => atype
    diheds::Vector{Pair{NTuple{4,Int16},Int16}} = Pair{NTuple{4,Int16},Int16}[] #dihedral is (i, j, k, l) => dtype
    improps::Vector{Pair{NTuple{4,Int16},Int16}} = Pair{NTuple{4,Int16},Int16}[] #improper is (i, j, k, l) => itype

    pair_coeffs::Vector{Tuple{Float64,Float64,Int16}} = Tuple{Float64,Float64,Int16}[]

    bond_coeffs::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]
    bond_map::Vector{Int16} = Int16[]

    angle_coeffs::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]
    angle_map::Vector{Int16} = Int16[]

    dihed_coeffs::Vector{NTuple{4,Float64}} = NTuple{4,Float64}[]
    dihed_map::Vector{Int16} = Int16[]

    improper_coeffs::Vector{Tuple{Float64,Int8,Int8}} = Tuple{Float64,Int8,Int8}[]
    improper_map::Vector{Int16} = Int16[]

    comment::String
end

"""
    read_lpg_data(f[; keywords...])

Read the datafile produced by LigParGen from `f`. `f` may be an I/O stream or a file name.

# Keywords
* `net_charge=nothing`: if set to a number, the charges will be tweaked so that the total
    charge equals to the specified value. If set to `nothing`, the charges will be tweaked
    so that the net charge is the nearest integer value
"""
function read_lpg_data(
    io::IO,
    ;
    net_charge::Union{Nothing,Real}=nothing,
)
    comment = readline(io)
    molstruct = Molecule(comment = comment[1] == '#' ? comment : "#" * comment)
    natom, nbond, nangle, ndihedral, nimproper,
    ntypes, nbtypes, natypes, ndtypes, nitypes = parse_header(io)

    resize!(molstruct.coords, natom)
    resize!(molstruct.types, natom)
    resize!(molstruct.typenames, natom)
    resize!(molstruct.masses, natom)
    resize!(molstruct.charges, natom)

    resize!(molstruct.bond_map, nbtypes)
    resize!(molstruct.angle_map, natypes)
    resize!(molstruct.dihed_map, ndtypes)
    resize!(molstruct.improper_map, nitypes)

    @assert startswith(readline(io), "Masses")

    read_masses!(molstruct, io, natom)

    readline(io)

    @assert startswith(readline(io), "Pair")

    read_pcoeff!(molstruct.pair_coeffs, io, ntypes)

    readline(io)

    @assert startswith(readline(io), "Bond")

    read_bcoeff!(molstruct.bond_coeffs, molstruct.bond_map, io, nbtypes)

    readline(io)

    @assert startswith(readline(io), "Angle")

    read_acoeff!(molstruct.angle_coeffs, molstruct.angle_map, io, natypes)
    readline(io)

    @assert startswith(readline(io), "Dihedral")

    read_dcoeff!(molstruct.dihed_coeffs, molstruct.dihed_map, io, ndtypes)

    readline(io)

    @assert startswith(readline(io), "Improper")

    read_icoeff!(molstruct.improper_coeffs, molstruct.improper_map, io, nitypes)

    readline(io)

    @assert startswith(readline(io), "Atoms")

    read_atoms!(molstruct.coords, molstruct.charges, molstruct.types, io, natom)

    readline(io)

    @assert startswith(readline(io), "Bonds")

    read_bonds!(molstruct.bonds, io, nbond, molstruct.bond_map)

    readline(io)

    @assert startswith(readline(io), "Angles")

    read_angles!(molstruct.angles, io, nangle, molstruct.angle_map)

    readline(io)

    @assert startswith(readline(io), "Dihedrals")

    read_dihed!(molstruct.diheds, io, ndihedral, molstruct.dihed_map)

    readline(io)

    @assert startswith(readline(io), "Impropers")

    read_impropers!(molstruct.improps, io, nimproper, molstruct.improper_map)

    if net_charge === nothing
        balance_charges!(molstruct, round(sum(molstruct.charges)); charge_diff_thresh=1.11e-3)
    else
        balance_charges!(molstruct, net_charge; charge_diff_thresh=1.11e-3)
    end

    refine_typenames!(molstruct)

    return molstruct
end

function read_lpg_data(fname::AbstractString; kw...)
    open(fname) do io
        read_lpg_data(io; kw...)
    end
end

function refine_typenames!(mol::Molecule)
    typenames = mol.typenames
    ct_types = (:CT1, :CT2, :CT3, :CT4)

    num_h_neighs = zero(mol.types)
    num_neighs = zero(mol.types)
    for bond in mol.bonds
        (i, k), = bond
        num_neighs[i] += true
        num_neighs[k] += true
        if typenames[i] == :H
            num_h_neighs[k] += true
        end
        if typenames[k] == :H
            num_h_neighs[i] += true
        end
    end
    for (k, nn) in enumerate(num_neighs)
        if typenames[k] == :C && nn == 4
            typenames[k] = :CT
        end
    end
    for (k, nh) in enumerate(num_h_neighs)
        if typenames[k] == :CT
            if nh > 0
                typenames[k] = ct_types[nh]
            else
                typenames[k] = :CT0
            end
        end
    end
    for bond in mol.bonds
        (i, k), = bond
        if typenames[i] == :H && typenames[k] in ct_types
            typenames[i] = :HC
        elseif typenames[k] == :H && typenames[i] in ct_types
            typenames[k] = :HC
        end
    end
    return mol
end
