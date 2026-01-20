const elems_by_mass = Dict(
    1 => :H,
    4 => :He,
    7 => :Li,
    9 => :Be,
    11 => :B,
    12 => :C,
    14 => :N,
    16 => :O,
    19 => :F,
    32 => :S,
    35 => :Cl,
    36 => :Cl,
)

Base.@kwdef struct LigParGenMolecule
    coords::Vector{NTuple{3,Float64}} = NTuple{3,Float64}[]
    types::Vector{Int16} = Int16[]
    charges::Vector{Float64} = Float64[]
    masses::Vector{Float64} = Float64[]
    bonds::Vector{NTuple{2,Int16}} = NTuple{2,Int16}[] #bond is (i, j)
    angles::Vector{NTuple{3,Int16}} = NTuple{3,Int16}[] #angle is (i, j, k)
    diheds::Vector{NTuple{4,Int16}} = NTuple{4,Int16}[] #dihedral is (i, j, k, l)
    improps::Vector{NTuple{4,Int16}} = NTuple{4,Int16}[] #improper is (i, j, k, l)

    pair_coeffs::Vector{Tuple{Float64,Float64,Int16}} = Tuple{Float64,Float64,Int16}[]
    bond_coeffs::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]
    angle_coeffs::Vector{NTuple{2,Float64}} = NTuple{2,Float64}[]
    dihed_coeffs::Vector{NTuple{4,Float64}} = NTuple{4,Float64}[]
    improper_coeffs::Vector{Tuple{Float64,Int8,Int8}} = Tuple{Float64,Int8,Int8}[]

    comment::String
    box_info::String
end

struct Molecule
    base::LigParGenMolecule
    compress_types::Bool
    compress_btypes::Bool
    compress_atypes::Bool
    compress_dtypes::Bool
    compress_itypes::Bool
    typenames::Vector{Symbol}
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
    natom, nbond, nangle, ndihedral, nimproper,
    ntypes, nbtypes, natypes, ndtypes, nitypes, box_info = parse_header(io)

    molstruct = LigParGenMolecule(
        ;
        comment = comment[1] == '#' ? comment : "#" * comment,
        box_info,
    )

    resize!(molstruct.coords, natom)
    resize!(molstruct.types, natom)
    resize!(molstruct.masses, natom)
    resize!(molstruct.charges, natom)

    @assert startswith(readline(io), "Masses")

    read_masses!(molstruct, io, natom)

    readline(io)

    @assert startswith(readline(io), "Pair")

    read_pcoeff!(molstruct.pair_coeffs, io, ntypes)

    readline(io)

    @assert startswith(readline(io), "Bond")

    read_bcoeff!(molstruct.bond_coeffs, io, nbtypes)

    readline(io)

    @assert startswith(readline(io), "Angle")

    read_acoeff!(molstruct.angle_coeffs, io, natypes)
    readline(io)

    @assert startswith(readline(io), "Dihedral")

    read_dcoeff!(molstruct.dihed_coeffs, io, ndtypes)

    readline(io)

    @assert startswith(readline(io), "Improper")

    read_icoeff!(molstruct.improper_coeffs, io, nitypes)

    readline(io)

    @assert startswith(readline(io), "Atoms")

    read_atoms!(molstruct.coords, molstruct.charges, molstruct.types, io, natom)

    readline(io)

    @assert startswith(readline(io), "Bonds")

    read_bonds!(molstruct.bonds, io, nbond)

    readline(io)

    @assert startswith(readline(io), "Angles")

    read_angles!(molstruct.angles, io, nangle)

    readline(io)

    @assert startswith(readline(io), "Dihedrals")

    read_dihed!(molstruct.diheds, io, ndihedral)

    readline(io)

    @assert startswith(readline(io), "Impropers")

    read_impropers!(molstruct.improps, io, nimproper)

    if net_charge === nothing
        balance_charges!(molstruct, round(sum(molstruct.charges)); charge_diff_thresh=1.11e-3)
    else
        balance_charges!(molstruct, net_charge; charge_diff_thresh=1.11e-3)
    end

    typenames = map(molstruct.masses) do m
        get(elems_by_mass, round(m), :X)
    end
    mol = Molecule(
        molstruct,
        compress_types,
        compress_btypes,
        compress_atypes,
        compress_dtypes,
        compress_itypes,
        typenames,
    )
    refine_typenames!(mol)

    return mol
end

function read_lpg_data(fname::AbstractString; kw...)
    open(fname) do io
        read_lpg_data(io; kw...)
    end
end

function refine_typenames!(molecule::Molecule)
    typenames = molecule.typenames

    mol = molecule.base
    ct_types = (:CT1, :CT2, :CT3, :CT4)

    num_h_neighs = zero(mol.types)
    num_neighs = zero(mol.types)
    for (i, k) in mol.bonds
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
    for (i, k) in mol.bonds
        if typenames[i] == :H && typenames[k] in ct_types
            typenames[i] = :HC
        elseif typenames[k] == :H && typenames[i] in ct_types
            typenames[k] = :HC
        end
    end
    return mol
end
