# Original contribution by Ivan M. Varaksin
# Refactoring by Vasily V. Pisarev

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

function fill_typenames!(mol::Molecule)
    typenames = mol.typenames
    ct_types = (:CT1, :CT2, :CT3, :CT4)

    for (k, amass) in enumerate(mol.masses)
        massnum = round(Int, amass)
        typenames[k] = elems_by_mass[massnum]
    end

    num_h_neighs = zero(mol.types)
    num_neighs = zero(mol.types)
    for bond in mol.bonds
        (i, k), = bond
        num_neighs[i] += true
        num_heighs[k] += true
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

"""
    switch_ff(mol::Molecule, ff_type::Symbol)

Switch LJ, bond, angle and dihedral parameters to the desired force field.

Known force fields: :opls_aa, :opls_aa_2020
"""
function switch_ff(mol_ref::Molecule, ff_type::Symbol)
    if !(ff_type in (:opls_aa, :opls_aa_2020))
        error("Unrecognized FF type, allowed types: :opls_aa, :opls_aa_2020")
    end

    mol = deepcopy(mol_ref)
    
    if ff_type == :opls_aa
        switch_ff_to_oplsaa!(mol)
    elseif ff_type == :opls_aa_2020
        switch_ff_to_opls2020!(mol)
    end
    return mol
end

"""
Скрипт для генерации OPLS потенциалов из LAMMPS .lmp файлов
Использование: julia opls_generator.jl <input.lmp> <potential_type>

Поддерживаемые типы потенциалов:
- opls-classic: OPLS-AA classic
- opls-2020: OPLS-AA/2020
- opls-classic-cm1a: OPLS-AA classic + CM1A
- opls-2020-cm1a: OPLS-AA/2020 + CM1A
"""

using Printf

# Структуры для хранения данных
struct AtomData
    id::Int
    mol_id::Int
    type_id::Int
    charge::Float64
    x::Float64
    y::Float64
    z::Float64
end

struct BondData
    id::Int
    type_id::Int
    atom1::Int
    atom2::Int
end

struct DihedralData
    id::Int
    type_id::Int
    atom1::Int
    atom2::Int
    atom3::Int
    atom4::Int
end

struct LammpsData
    masses::Dict{Int, Float64}
    atoms::Vector{AtomData}
    bonds::Vector{BondData}
    dihedrals::Vector{DihedralData}
    pair_coeffs::Dict{Int, Tuple{Float64, Float64}}
    bond_coeffs::Dict{Int, Tuple{Float64, Float64}}
    dihedral_coeffs::Dict{Int, Tuple{Float64, Float64, Float64, Float64}}
end

"""
Парсинг .lmp файла
"""
function parse_lmp_file(filename::String)::LammpsData
    masses = Dict{Int, Float64}()
    atoms = Vector{AtomData}()
    bonds = Vector{BondData}()
    dihedrals = Vector{DihedralData}()
    pair_coeffs = Dict{Int, Tuple{Float64, Float64}}()
    bond_coeffs = Dict{Int, Tuple{Float64, Float64}}()
    dihedral_coeffs = Dict{Int, Tuple{Float64, Float64, Float64, Float64}}()

    current_section = ""

    open(filename, "r") do file
        for line in eachline(file)
            line = strip(line)

            # Пропускаем пустые строки и комментарии
            if isempty(line) || startswith(line, "#")
                continue
            end

            # Определяем секции
            if line == "Masses"
                current_section = "masses"
                continue
            elseif line == "Pair Coeffs"
                current_section = "pair_coeffs"
                continue
            elseif line == "Bond Coeffs"
                current_section = "bond_coeffs"
                continue
            elseif line == "Dihedral Coeffs"
                current_section = "dihedral_coeffs"
                continue
            elseif line == "Atoms"
                current_section = "atoms"
                continue
            elseif line == "Bonds"
                current_section = "bonds"
                continue
            elseif line == "Dihedrals"
                current_section = "dihedrals"
                continue
            elseif line == "Angles"
                current_section = "angles"
                continue
            end

            # Парсим данные в зависимости от секции
            parts = split(line)

            if current_section == "masses" && length(parts) >= 2
                type_id = parse(Int, parts[1])
                mass = parse(Float64, parts[2])
                masses[type_id] = mass

            elseif current_section == "pair_coeffs" && length(parts) >= 3
                type_id = parse(Int, parts[1])
                epsilon = parse(Float64, parts[2])
                sigma = parse(Float64, parts[3])
                pair_coeffs[type_id] = (epsilon, sigma)

            elseif current_section == "bond_coeffs" && length(parts) >= 3
                type_id = parse(Int, parts[1])
                k = parse(Float64, parts[2])
                r0 = parse(Float64, parts[3])
                bond_coeffs[type_id] = (k, r0)

            elseif current_section == "dihedral_coeffs" && length(parts) >= 5
                type_id = parse(Int, parts[1])
                k1 = parse(Float64, parts[2])
                k2 = parse(Float64, parts[3])
                k3 = parse(Float64, parts[4])
                k4 = parse(Float64, parts[5])
                dihedral_coeffs[type_id] = (k1, k2, k3, k4)

            elseif current_section == "atoms" && length(parts) >= 7
                id = parse(Int, parts[1])
                mol_id = parse(Int, parts[2])
                type_id = parse(Int, parts[3])
                charge = parse(Float64, parts[4])
                x = parse(Float64, parts[5])
                y = parse(Float64, parts[6])
                z = parse(Float64, parts[7])
                push!(atoms, AtomData(id, mol_id, type_id, charge, x, y, z))

            elseif current_section == "bonds" && length(parts) >= 4
                id = parse(Int, parts[1])
                type_id = parse(Int, parts[2])
                atom1 = parse(Int, parts[3])
                atom2 = parse(Int, parts[4])
                push!(bonds, BondData(id, type_id, atom1, atom2))

            elseif current_section == "angles" && length(parts) >= 5
                # Из углов извлекаем связи: угол A-B-C содержит связи A-B и B-C
                atom1 = parse(Int, parts[3])
                atom2 = parse(Int, parts[4])
                atom3 = parse(Int, parts[5])
                # Добавляем связи A-B и B-C (используем фиктивные id и type)
                push!(bonds, BondData(length(bonds) + 1, 1, atom1, atom2))
                push!(bonds, BondData(length(bonds) + 1, 1, atom2, atom3))

            elseif current_section == "dihedrals" && length(parts) >= 6
                id = parse(Int, parts[1])
                type_id = parse(Int, parts[2])
                atom1 = parse(Int, parts[3])
                atom2 = parse(Int, parts[4])
                atom3 = parse(Int, parts[5])
                atom4 = parse(Int, parts[6])
                push!(dihedrals, DihedralData(id, type_id, atom1, atom2, atom3, atom4))
            end
        end
    end

    # Удаляем дублирующиеся связи (из углов могут получиться повторы)
    unique_bonds = Vector{BondData}()
    seen_bonds = Set{Tuple{Int, Int}}()

    for bond in bonds
        # Нормализуем порядок атомов (меньший id первым)
        atom_pair = bond.atom1 < bond.atom2 ? (bond.atom1, bond.atom2) : (bond.atom2, bond.atom1)
        if !(atom_pair in seen_bonds)
            push!(seen_bonds, atom_pair)
            push!(unique_bonds, bond)
        end
    end

    return LammpsData(masses, atoms, unique_bonds, dihedrals, pair_coeffs, bond_coeffs, dihedral_coeffs)
end

"""
Определение типов атомов (углерод/водород) по массе
"""
function identify_atom_types(data::LammpsData)
    carbon_types = Set{Int}()
    hydrogen_types = Set{Int}()

    for (type_id, mass) in data.masses
        if round(mass) == 12
            push!(carbon_types, type_id)
        elseif round(mass) == 1
            push!(hydrogen_types, type_id)
        end
    end

    return carbon_types, hydrogen_types
end

"""
Построение графа связей для анализа структуры
"""
function build_bond_graph(data::LammpsData)
    # Создаем словарь: атом -> список связанных атомов
    bond_graph = Dict{Int, Vector{Int}}()

    # Инициализируем пустые списки для всех атомов
    for atom in data.atoms
        bond_graph[atom.id] = Int[]
    end

    # Добавляем связи
    println("Отладка связей:")
    println("  Всего связей: $(length(data.bonds))")
    if length(data.bonds) > 0
        println("  Первые 5 связей:")
        for i in 1:min(5, length(data.bonds))
            bond = data.bonds[i]
            println("    Связь $i: атом $(bond.atom1) <-> атом $(bond.atom2)")
        end
    end

    for bond in data.bonds
        push!(bond_graph[bond.atom1], bond.atom2)
        push!(bond_graph[bond.atom2], bond.atom1)
    end

    # Проверяем несколько атомов
    println("  Связи для первых 3 атомов:")
    for i in 1:min(3, length(data.atoms))
        atom = data.atoms[i]
        neighbors = get(bond_graph, atom.id, Int[])
        println("    Атом $(atom.id) (тип $(atom.type_id)): связан с $(neighbors)")
    end

    return bond_graph
end

"""
Классификация углеродов на терминальные и средние
"""
function classify_carbons(data::LammpsData, carbon_types::Set{Int}, hydrogen_types::Set{Int})
    bond_graph = build_bond_graph(data)
    terminal_carbon_types = Set{Int}()
    midchain_carbon_types = Set{Int}()

    # Создаем словарь: id атома -> тип атома
    atom_type_map = Dict{Int, Int}()
    for atom in data.atoms
        atom_type_map[atom.id] = atom.type_id
    end

    # Подсчитываем связи для каждого типа углерода
    carbon_hydrogen_bonds = Dict{Int, Int}()

    # Анализируем каждый углерод
    for atom in data.atoms
        if atom.type_id in carbon_types
            # Считаем количество связей с водородами для этого атома
            hydrogen_bonds = 0
            if haskey(bond_graph, atom.id)
                for neighbor_id in bond_graph[atom.id]
                    if haskey(atom_type_map, neighbor_id)
                        neighbor_type = atom_type_map[neighbor_id]
                        if neighbor_type in hydrogen_types
                            hydrogen_bonds += 1
                        end
                    end
                end
            end

            # Сохраняем максимальное количество связей для этого типа углерода
            if !haskey(carbon_hydrogen_bonds, atom.type_id) || hydrogen_bonds > carbon_hydrogen_bonds[atom.type_id]
                carbon_hydrogen_bonds[atom.type_id] = hydrogen_bonds
            end
        end
    end

    # Классифицируем типы углеродов по количеству связей с водородами
    for (carbon_type, hydrogen_bonds) in carbon_hydrogen_bonds
        if hydrogen_bonds == 3
            push!(terminal_carbon_types, carbon_type)
        elseif hydrogen_bonds == 2
            push!(midchain_carbon_types, carbon_type)
        end
    end

    println("Отладка классификации углеродов:")
    println("  Связи углерод-водород: $carbon_hydrogen_bonds")
    println("  Терминальные углероды: $terminal_carbon_types")
    println("  Средние углероды: $midchain_carbon_types")

    return terminal_carbon_types, midchain_carbon_types
end

"""
Поиск диэдральных углов C-C-C-C
"""
function find_cccc_dihedrals(data::LammpsData, carbon_types::Set{Int})
    cccc_dihedral_types = Set{Int}()

    # Создаем словарь: id атома -> тип атома
    atom_type_map = Dict{Int, Int}()
    for atom in data.atoms
        atom_type_map[atom.id] = atom.type_id
    end

    # Проверяем каждый диэдральный угол
    for dihedral in data.dihedrals
        # Проверяем, что все четыре атома - углероды
        if (atom_type_map[dihedral.atom1] in carbon_types &&
            atom_type_map[dihedral.atom2] in carbon_types &&
            atom_type_map[dihedral.atom3] in carbon_types &&
            atom_type_map[dihedral.atom4] in carbon_types)
            push!(cccc_dihedral_types, dihedral.type_id)
        end
    end

    return cccc_dihedral_types
end

"""
Генерация команд для OPLS-AA classic
"""
function generate_opls_classic(data::LammpsData, carbon_types::Set{Int}, hydrogen_types::Set{Int},
                              terminal_carbon_types::Set{Int}, midchain_carbon_types::Set{Int})
    commands = String[]

    # Заряды
    push!(commands, "# Charges")
    for (type_id, _) in data.masses
        if type_id in hydrogen_types
            push!(commands, "set type $type_id charge 0.06")
        elseif type_id in terminal_carbon_types
            push!(commands, "set type $type_id charge -0.18")
        elseif type_id in midchain_carbon_types
            push!(commands, "set type $type_id charge -0.12")
        end
    end

    push!(commands, "")
    push!(commands, "# Pair coefficients (original FF)")
    for (type_id, (epsilon, sigma)) in data.pair_coeffs
        push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, epsilon, sigma))
    end

    push!(commands, "")
    push!(commands, "# Bond coefficients (original FF)")
    for (type_id, (k, r0)) in data.bond_coeffs
        push!(commands, @sprintf("bond_coeff %2d %.4f %.4f", type_id, k, r0))
    end

    push!(commands, "")
    push!(commands, "# Dihedral coefficients (original FF)")
    for (type_id, (k1, k2, k3, k4)) in data.dihedral_coeffs
        push!(commands, @sprintf("dihedral_coeff %2d %.3f %.3f %.3f %.3f", type_id, k1, k2, k3, k4))
    end

    return commands
end

"""
Генерация команд для OPLS-AA/2020
"""
function generate_opls_2020(data::LammpsData, carbon_types::Set{Int}, hydrogen_types::Set{Int},
                           terminal_carbon_types::Set{Int}, midchain_carbon_types::Set{Int},
                           cccc_dihedral_types::Set{Int})
    commands = String[]

    # Заряды
    push!(commands, "# Charges")
    for (type_id, _) in data.masses
        if type_id in hydrogen_types
            push!(commands, "set type $type_id charge 0.06")
        elseif type_id in terminal_carbon_types
            push!(commands, "set type $type_id charge -0.18")
        elseif type_id in midchain_carbon_types
            push!(commands, "set type $type_id charge -0.12")
        end
    end

    push!(commands, "")
    push!(commands, "# Pair coefficients (OPLS-AA/2020)")
    for (type_id, _) in data.masses
        if type_id in hydrogen_types
            push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, 0.026, 2.480))
        elseif type_id in terminal_carbon_types
            push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, 0.066, 3.550))
        elseif type_id in midchain_carbon_types
            push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, 0.066, 3.510))
        end
    end

    push!(commands, "")
    push!(commands, "# Bond coefficients (original FF)")
    for (type_id, (k, r0)) in data.bond_coeffs
        push!(commands, @sprintf("bond_coeff %2d %.4f %.4f", type_id, k, r0))
    end

    push!(commands, "")
    push!(commands, "# Dihedral coefficients")
    for (type_id, _) in data.dihedral_coeffs
        if type_id in cccc_dihedral_types
            push!(commands, @sprintf("dihedral_coeff %2d %.2f %.2f %.2f %.1f", type_id, 0.85, -0.20, 0.20, 0.0))
        else
            # Оригинальные коэффициенты для не C-C-C-C диэдралов
            k1, k2, k3, k4 = data.dihedral_coeffs[type_id]
            push!(commands, @sprintf("dihedral_coeff %2d %.3f %.3f %.3f %.3f", type_id, k1, k2, k3, k4))
        end
    end

    return commands
end

"""
Генерация команд для OPLS-AA classic + CM1A
"""
function generate_opls_classic_cm1a(data::LammpsData, carbon_types::Set{Int}, hydrogen_types::Set{Int},
                                   terminal_carbon_types::Set{Int}, midchain_carbon_types::Set{Int})
    commands = String[]

    # Заряды НЕ устанавливаем - используем оригинальные из LigParGen
    push!(commands, "# Charges from LigParGen (CM1A) - already set in data file")
    push!(commands, "")

    push!(commands, "# Pair coefficients (original FF)")
    for (type_id, (epsilon, sigma)) in data.pair_coeffs
        push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, epsilon, sigma))
    end

    push!(commands, "")
    push!(commands, "# Bond coefficients (original FF)")
    for (type_id, (k, r0)) in data.bond_coeffs
        push!(commands, @sprintf("bond_coeff %2d %.4f %.4f", type_id, k, r0))
    end

    push!(commands, "")
    push!(commands, "# Dihedral coefficients (original FF)")
    for (type_id, (k1, k2, k3, k4)) in data.dihedral_coeffs
        push!(commands, @sprintf("dihedral_coeff %2d %.3f %.3f %.3f %.3f", type_id, k1, k2, k3, k4))
    end

    return commands
end

"""
Генерация команд для OPLS-AA/2020 + CM1A
"""
function generate_opls_2020_cm1a(data::LammpsData, carbon_types::Set{Int}, hydrogen_types::Set{Int},
                                terminal_carbon_types::Set{Int}, midchain_carbon_types::Set{Int},
                                cccc_dihedral_types::Set{Int})
    commands = String[]

    # Заряды НЕ устанавливаем - используем оригинальные из LigParGen
    push!(commands, "# Charges from LigParGen (CM1A) - already set in data file")
    push!(commands, "")

    push!(commands, "# Pair coefficients (OPLS-AA/2020)")
    for (type_id, _) in data.masses
        if type_id in hydrogen_types
            push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, 0.026, 2.480))
        elseif type_id in terminal_carbon_types
            push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, 0.066, 3.550))
        elseif type_id in midchain_carbon_types
            push!(commands, @sprintf("pair_coeff %2d %2d %.3f %.7f", type_id, type_id, 0.066, 3.510))
        end
    end

    push!(commands, "")
    push!(commands, "# Bond coefficients (original FF)")
    for (type_id, (k, r0)) in data.bond_coeffs
        push!(commands, @sprintf("bond_coeff %2d %.4f %.4f", type_id, k, r0))
    end

    push!(commands, "")
    push!(commands, "# Dihedral coefficients")
    for (type_id, _) in data.dihedral_coeffs
        if type_id in cccc_dihedral_types
            push!(commands, @sprintf("dihedral_coeff %2d %.2f %.2f %.2f %.1f", type_id, 0.85, -0.20, 0.20, 0.0))
        else
            # Оригинальные коэффициенты для не C-C-C-C диэдралов
            k1, k2, k3, k4 = data.dihedral_coeffs[type_id]
            push!(commands, @sprintf("dihedral_coeff %2d %.3f %.3f %.3f %.3f", type_id, k1, k2, k3, k4))
        end
    end

    return commands
end

"""
Основная функция
"""
function main()
    if length(ARGS) != 2
        println("Использование: julia opls_generator.jl <input.lmp> <potential_type>")
        println("Поддерживаемые типы потенциалов:")
        println("  opls-classic      - OPLS-AA classic")
        println("  opls-2020         - OPLS-AA/2020")
        println("  opls-classic-cm1a - OPLS-AA classic + CM1A")
        println("  opls-2020-cm1a    - OPLS-AA/2020 + CM1A")
        exit(1)
    end

    input_file = ARGS[1]
    potential_type = ARGS[2]

    if !isfile(input_file)
        println("Ошибка: файл $input_file не найден")
        exit(1)
    end

    # Парсим входной файл
    println("Парсинг файла $input_file...")
    data = parse_lmp_file(input_file)

    # Анализируем структуру
    println("Анализ структуры молекулы...")
    carbon_types, hydrogen_types = identify_atom_types(data)
    terminal_carbon_types, midchain_carbon_types = classify_carbons(data, carbon_types, hydrogen_types)
    cccc_dihedral_types = find_cccc_dihedrals(data, carbon_types)

    println("Найдено:")
    println("  Типы углеродов: $(length(carbon_types))")
    println("  Типы водородов: $(length(hydrogen_types))")
    println("  Терминальные углероды: $(length(terminal_carbon_types))")
    println("  Средние углероды: $(length(midchain_carbon_types))")
    println("  C-C-C-C диэдралы: $(length(cccc_dihedral_types))")

    # Генерируем команды
    commands = String[]
    if potential_type == "opls-classic"
        commands = generate_opls_classic(data, carbon_types, hydrogen_types,
                                       terminal_carbon_types, midchain_carbon_types)
    elseif potential_type == "opls-2020"
        commands = generate_opls_2020(data, carbon_types, hydrogen_types,
                                    terminal_carbon_types, midchain_carbon_types, cccc_dihedral_types)
    elseif potential_type == "opls-classic-cm1a"
        commands = generate_opls_classic_cm1a(data, carbon_types, hydrogen_types,
                                            terminal_carbon_types, midchain_carbon_types)
    elseif potential_type == "opls-2020-cm1a"
        commands = generate_opls_2020_cm1a(data, carbon_types, hydrogen_types,
                                         terminal_carbon_types, midchain_carbon_types, cccc_dihedral_types)
    else
        println("Ошибка: неизвестный тип потенциала $potential_type")
        exit(1)
    end

    # Выводим результат
    output_file = replace(input_file, ".lmp" => "_$(potential_type).in")
    open(output_file, "w") do file
        for command in commands
            println(file, command)
        end
    end

    println("Команды сохранены в файл: $output_file")
end

# Запуск основной функции
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
