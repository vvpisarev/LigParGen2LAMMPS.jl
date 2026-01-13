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
* `compress_types::Bool=false`: mark particles with the same LJ parameters as the same type
* `compress_btypes::Bool=true`: mark bonds with the same parameters as the same type
* `compress_atypes::Bool=true`: mark angles with the same parameters as the same type
* `compress_dtypes::Bool=true`: mark dihedrals with the same parameters as the same type
* `compress_itypes::Bool=true`: mark impropers with the same parameters as the same type
* `net_charge=nothing`: if set to a number, the charges will be tweaked so that the total
    charge equals to the specified value. If set to `nothing`, the charges will be tweaked
    so that the net charge is the nearest integer value
"""
function read_lpg_data(
    io::IO,
    ;
    compress_types::Bool=false,
    compress_btypes::Bool=true,
    compress_atypes::Bool=true,
    compress_dtypes::Bool=true,
    compress_itypes::Bool=true,
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

    read_masses!(molstruct.masses, io, natom)

    readline(io)

    @assert startswith(readline(io), "Pair")

    read_pcoeff!(molstruct.pair_coeffs, io, ntypes)

    readline(io)

    @assert startswith(readline(io), "Bond")

    read_bcoeff!(molstruct.bond_coeffs, molstruct.bond_map, io, nbtypes, compress_btypes)

    readline(io)

    @assert startswith(readline(io), "Angle")

    read_acoeff!(molstruct.angle_coeffs, molstruct.angle_map, io, natypes, compress_atypes)

    readline(io)

    @assert startswith(readline(io), "Dihedral")

    read_dcoeff!(molstruct.dihed_coeffs, molstruct.dihed_map, io, ndtypes, compress_dtypes)

    readline(io)

    @assert startswith(readline(io), "Improper")

    read_icoeff!(molstruct.improper_coeffs, molstruct.improper_map, io, nitypes, compress_itypes)

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

    fill_typenames!(molstruct)

    if compress_types
        compress_types!(molstruct)
    end

    return molstruct
end

function read_lpg_data(fname::AbstractString; kw...)
    open(fname) do io
        read_lpg_data(io; kw...)
    end
end
