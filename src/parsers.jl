function parse_header(io)
    natom = nbond = nangle = ndihedral = nimproper = 0
    ntypes = nbtypes = natypes = ndtypes = nitypes = 0

    boxinfo_buf = UInt8[]
    while !isuppercase(peek(io, Char))
        datastr = readline(io) |> strip
        isempty(datastr) && continue
        if endswith(datastr, "atoms")
            natom = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "bonds")
            nbond = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "angles")
            nangle = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "dihedrals")
            ndihedral = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "impropers")
            nimproper = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "atom types")
            ntypes = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "bond types")
            nbtypes = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "angle types")
            natypes = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "dihedral types")
            ndtypes = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "improper types")
            nitypes = parse(Int, rstrip(!isdigit, datastr))
        elseif endswith(datastr, "lo")
            append!(boxinfo_buf, codeunits(datastr), ('\n',))
        else
            @warn "Unused header line: $datastr"
        end
    end
    natom, nbond, nangle, ndihedral, nimproper,
    ntypes, nbtypes, natypes, ndtypes, nitypes, String(boxinfo_buf)
end

function read_masses!(mol::LigParGenMolecule, io::IO, natoms::Integer)
    (; masses) = mol
    readline(io)

    for _ in 1:natoms
        ln = readline(io)
        if isempty(ln)
            error("Not all masses are given")
        end
        is, ms = split(ln)
        i = parse(Int16, is)
        m = parse(Float64, ms)
        masses[i] = m
    end
    return mol
end

function read_pcoeff!(pair_coeff, io::IO, ntypes::Integer)
    readline(io)

    for ia in 1:ntypes
        ln = readline(io)
        if isempty(ln)
            error("Not all pair coeffs are given")
        end
        is, εs, σs = split(ln)
        itype = parse(Int16, is)
        ε = parse(Float64, εs)
        σ = parse(Float64, σs)
        push!(pair_coeff, (ε, σ, itype))
    end
    return pair_coeff
end

function read_bcoeff!(bond_coeff, io::IO, nbonds::Integer)
    readline(io)

    for ib in 1:nbonds
        ln = readline(io)
        if isempty(ln)
            error("Not all bond coeffs are given")
        end
        is, ks, rs = split(ln)
        i = parse(Int16, is)
        k = parse(Float64, ks)
        r0 = parse(Float64, rs)
        push!(bond_coeff, (k, r0))
    end
    return bond_coeff
end

function read_acoeff!(angle_coeff, io::IO, nangles::Integer)
    readline(io)

    for ia in 1:nangles
        ln = readline(io)
        if isempty(ln)
            error("Not all angle coeffs are given")
        end
        is, ks, θs = split(ln)
        i = parse(Int16, is)
        k = parse(Float64, ks)
        θ0 = parse(Float64, θs)
        push!(angle_coeff, (k, θ0))
    end
    return angle_coeff
end

function read_dcoeff!(dihed_coeff, io::IO, ndihed::Integer)
    readline(io)

    for id in 1:ndihed
        ln = readline(io)
        if isempty(ln)
            error("Not all dihedral coeffs are given")
        end
        is, c1s, c2s, c3s, c4s = split(ln)
        i = parse(Int16, is)
        c1, c2, c3, c4 = parse.(Float64, (c1s, c2s, c3s, c4s))

        push!(dihed_coeff, (c1, c2, c3, c4))
    end
    return dihed_coeff
end

function read_icoeff!(improper_coeff, io::IO, nimproper::Integer)
    readline(io)

    for ii in 1:nimproper
        ln = readline(io)
        if isempty(ln)
            error("Not all improper coeffs are given")
        end
        is, ks, ds, ns = split(ln)
        i = parse(Int16, is)
        k = parse(Float64, ks)
        d = parse(Int8, ds)
        n = parse(Int8, ns)

        push!(improper_coeff, (k, d, n))
    end
    return improper_coeff
end

function read_atoms!(coord, charge, type, io::IO, natoms::Integer)
    readline(io)

    for _ in 1:natoms
        ln = readline(io)
        if isempty(ln)
            error("Missing atoms")
        end
        is, ms, ts, qs, xs, ys, zs = split(ln)
        i, t = parse.(Int16, (is, ts))
        q, x, y, z = parse.(Float64, (qs, xs, ys, zs))
        charge[i] = q
        coord[i] = (x, y, z)
        type[i] = t
    end

    return coord, charge
end

function read_bonds!(bonds, io::IO, nbonds::Integer)
    readline(io)

    for _ in 1:nbonds
        ln = readline(io)
        if isempty(ln)
            error("Missing bonds")
        end
        is, ts, a1s, a2s = split(ln)
        _, type, a1, a2 = parse.(Int16, (is, ts, a1s, a2s))
        push!(bonds, (a1, a2))
    end
    return bonds
end

function read_angles!(angles, io::IO, nangles::Integer)
    readline(io)

    for _ in 1:nangles
        ln = readline(io)
        if isempty(ln)
            error("Missing angles")
        end
        is, ts, a1s, a2s, a3s = split(ln)
        i, type, a1, a2, a3 = parse.(Int16, (is, ts, a1s, a2s, a3s))
        push!(angles, (a1, a2, a3))
    end
    return angles
end

function read_dihed!(dihed, io::IO, ndihed::Integer)
    readline(io)

    for _ in 1:ndihed
        ln = readline(io)
        if isempty(ln)
            error("Missing dihedrals")
        end
        is, ts, a1s, a2s, a3s, a4s = split(ln)
        i, type, a1, a2, a3, a4 = parse.(Int16, (is, ts, a1s, a2s, a3s, a4s))
        push!(dihed, (a1, a2, a3, a4))
    end
    return dihed
end

function read_impropers!(improper, io::IO, nimproper::Integer)
    readline(io)

    for _ in 1:nimproper
        ln = readline(io)
        if isempty(ln)
            error("Missing impropers")
        end
        is, ts, a1s, a2s, a3s, a4s = split(ln)
        i, type, a1, a2, a3, a4 = parse.(Int16, (is, ts, a1s, a2s, a3s, a4s))
        push!(improper, (a1, a2, a3, a4))
    end
    return improper
end
