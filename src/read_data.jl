export read_lpg_data

include("parsers.jl")

Base.@kwdef struct Molecule
    coords::Vector{NTuple{3,Float64}} = NTuple{3,Float64}[]
    types::Vector{Int16} = Int16[]
    charges::Vector{Float64} = Float64[]
    masses::Vector{Float64} = Float64[]
    bonds::Vector{Pair{NTuple{2,Int16},Int16}} = Pair{NTuple{2,Int16},Int16}[] #bond is (i, j) => btype
    angles::Vector{Pair{NTuple{3,Int16},Int16}} = Pair{NTuple{3,Int16},Int16}[] #angle is (i, j, k) => atype
    diheds::Vector{Pair{NTuple{4,Int16},Int16}} = Pair{NTuple{4,Int16},Int16}[] #dihedral is (i, j, k, l) => dtype
    improps::Vector{Pair{NTuple{4,Int16},Int16}} = Pair{NTuple{4,Int16},Int16}[] #improper is (i, j, k, l) => itype

    pair_coeffs::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]

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
    read_lpg_data(f)

Read the datafile produced by LigParGen from `f`. `f` may 
be an I/O stream or a file name.
"""
function read_lpg_data(io::IO)
    comment = readline(io)
    molstruct = Molecule(comment = comment[1] == '#' ? comment : "#" * comment)
    natom, nbond, nangle, ndihedral, nimproper,
    ntypes, nbtypes, natypes, ndtypes, nitypes = parse_header(io)

    resize!(molstruct.coords, natom)
    resize!(molstruct.types, natom)
    resize!(molstruct.masses, natom)
    resize!(molstruct.charges, natom)

    resize!(molstruct.bond_map, nbtypes)
    resize!(molstruct.angle_map, natypes)
    resize!(molstruct.dihed_map, ndtypes)
    resize!(molstruct.improper_map, nitypes)

    @assert startswith(readline(io), "Masses")

    read_masses!(molstruct.masses, io, natom)

    readline(io)

    @assert startswith(readline(io), "Pair")

    read_pcoeff!(molstruct.pair_coeffs, molstruct.types, io, ntypes)

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

    read_atoms!(molstruct.coords, molstruct.charges, io, natom)

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

    return molstruct
end

read_lpg_data(fname::AbstractString) = open(fname) do io read_lpg_data(io) end