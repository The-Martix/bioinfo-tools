# This script/library gets the PDB structure of a given pdb file and set structural classifications measuring distances of CAMs within a given threshold

from Bio.PDB import PDBParser, MMCIFIO, MMCIFParser
from Bio.PDB import Superimposer
from Bio.PDB import PDBIO, StructureBuilder
from Bio.PDB.ResidueDepth import get_surface as GetRDSurface
from Bio.PDB.ResidueDepth import min_dist
import copy
import aminotools
import numpy as np
from scipy.spatial import cKDTree
import os

residues_names = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


############################################################################### PROCESAMIENTO DE PDB FILES ############################################################################################

# Remueve elementos vacios de una lista
def remove_gaps(lista):
    new_list = []
    for elem in lista: 
        if elem != "" and elem != None: new_list.append(elem.replace("\n", "").strip())
    return new_list

# Esto es para generar espacios vacios en un string, lo uso para escribir bien la linea del pdb
def generate_spaces(amount):
    return " "*amount

# Obtengo las lineas de un .pdb
def open_pdb(pdb_file):
    with open(pdb_file, "r") as file:
        lines = file.readlines()
    return lines

# Obtener los campos de una linea de ATOM o HETATOM de PDB
def parse_pdb_line(line):
    return {
        "record": line[0:6].strip(),
        "atom_number": int(line[6:11]),
        "atom_name": line[12:16].strip(),
        "residue_name": line[17:20].strip(),
        "chain_id": line[21],
        "residue_number": int(line[22:26]),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "occupancy": float(line[54:60]),
        "temp_factor": float(line[60:66]),
        "element": line[76:78].strip()
    }

# Volver a reconstruir la linea del PDB parseada
def unparse_pdb_line(atom_dict):
    return (
        f"{atom_dict['record']:<6}"              # Record name (ATOM/HETATM)
        f"{int(atom_dict['serial']):>5} "         # Atom serial number
        f"{atom_dict['name']:>4}"                 # Atom name
        f"{atom_dict['altLoc'] or ' ':1}"         # Alternate location indicator
        f"{atom_dict['resName']:<3} "             # Residue name
        f"{atom_dict['chainID']:>1}"              # Chain ID
        f"{int(atom_dict['resSeq']):>4}"          # Residue sequence number
        f"    "                                   # 4 spaces (insertion code unused)
        f"{float(atom_dict['x']):>8.3f}"          # X
        f"{float(atom_dict['y']):>8.3f}"          # Y
        f"{float(atom_dict['z']):>8.3f}"          # Z
        f"{float(atom_dict['occupancy']):>6.2f}"  # Occupancy
        f"{float(atom_dict['tempFactor']):>6.2f}" # Temp factor
        f"          {atom_dict['element']:>2}"    # Element
        f"\n"
    )

# Escribir archivo pdb
def write_pdb(pdb_file, new_lines, show_msg=True):
    file = open(pdb_file, 'w')
    for line in new_lines: file.write(line)
    file.close()
    if show_msg: print(f"Se ha escrito el archivo {pdb_file} con exito.")

# Esta funcion te deja cambiar el identifier de una linea de atomo (por ej ATOM->HETATM, etc)
def update_identifier(lines, old_identifier, new_identifier):
    new_lines = []
    for line in lines:
        if line.startswith(old_identifier):
            line = line.replace(old_identifier, new_identifier)
        new_lines.append(line)
    return new_lines

# Obtengo las lineas de atomos del ligando del Docking
def get_hetatoms_lines(lines, hetatoms:list[str]):
    atom_lines = []
    for line in lines:
        if line.startswith("HETATOM") or line.startswith("HETATM"):
            splits = remove_gaps(line.split(" "))
            if len(splits) < 4: continue  # Evitar lineas mal form
            hetatom_id = splits[3]
            if hetatom_id in hetatoms:
                atom_lines.append(line)
    return atom_lines

# Esta funcion toma un PDB en el que estan tus ligandos de interes, una lista con los nombres de los ligandos y devuelve un pdb con solo los ligandos de interes
def write_ligands_pdb(input_pdb, ligands: list[str], output_pdb):
    lines = open_pdb(input_pdb)
    ligand_lines = get_hetatoms_lines(lines, hetatoms=ligands)
    write_pdb(output_pdb, ligand_lines)

# Renombrar ligando
def rename_hetatom(pdb_lines, hetatom_dict):
    new_lines = []
    for line in pdb_lines:
        if line.startswith("HETATOM") or line.startswith("HETATM"):
            for key in list(hetatom_dict.keys()):
                if key in line:
                    line = line.replace(key, hetatom_dict[key])
                    break
        new_lines.append(line)
    return new_lines

# Modificar la linea del ligando para que sea compatible con la del pdb de la proteina
def update_ligand_line(ligand_line, last_reference_line):
    last_line_splits = remove_gaps(last_reference_line.split(" "))
    last_atom_index = int(last_line_splits[1])
    last_atom_chain = str(last_line_splits[4])
    last_resd_index = int(last_line_splits[5])

    splits = remove_gaps(ligand_line.split(" "))
    atom_index = str(last_atom_index + int(splits[1]))
    resd_index = str(last_resd_index + int(splits[4]))
    new_line = f"{splits[0]} {atom_index}{generate_spaces(6-len(atom_index))}{splits[2]}{generate_spaces(4-len(splits[2]))}{splits[3]}{generate_spaces(4-len(splits[3]))}{last_atom_chain} {resd_index}{generate_spaces(9-len(resd_index))}{splits[5]}  {splits[6]}  {splits[7]}  {splits[8]} {splits[9]}{generate_spaces(16-len(splits[9]))}{splits[10]}\n"
    return new_line

    
# Meto las lineas del ligando en el pdb de la proteina
def merge_ligand_to_pdb(reference_pdb, ligand_pdb, hetatoms):
    atom_lines = get_hetatoms_lines(open_pdb(ligand_pdb), hetatoms=hetatoms)
    reference_lines = open_pdb(reference_pdb)
    new_lines = []
    for line in reference_lines[:-1]: new_lines.append(line)
    for lig_line in atom_lines:
        new_line = update_ligand_line(ligand_line=lig_line, last_reference_line=line)
        new_lines.append(new_line)
    new_lines.append(reference_lines[-1])
    return new_lines

# Borrar HETATOMS del pdb
def filter_hetatoms(lines, accepted_hetatoms=[]):
    new_lines = []    
    for line in lines:
        if line.startswith("HETATM") or line.startswith("ANISOU"):
            splits = line.split(" ")
            for i in range(splits.count("")): splits.remove("")
            hetatom_id = splits[3]
            if hetatom_id in accepted_hetatoms: new_lines.append(line)
        elif line.startswith("HET"):
            splits = line.split(" ")
            for i in range(splits.count("")): splits.remove("")
            hetatom_id = splits[1]
            if hetatom_id in accepted_hetatoms: new_lines.append(line)
        else: new_lines.append(line)
    return new_lines

# Borrar aguas del PDB
def remove_waters(lines):
    new_lines = []
    for line in lines:
        if "HETATOM" in line or "HETATM" in line:
            if "HOH" in line: continue
        new_lines.append(line)
    return new_lines

# Borrar cadenas del pdb
def filter_chains(lines, accepted_chains=[]):
    new_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
            splits = line.split(" ")
            for i in range(splits.count("")): splits.remove("")
            chain_id = splits[4]
            if chain_id in accepted_chains: new_lines.append(line)
        else: new_lines.append(line)
    return new_lines

# Reindexar las posiciones de cada cadena para que cada una empiece desde el 1 (o el start_index que se pase como input)
def reindex_chains(lines, start_index=1):
    new_lines, chain_map = [], {}
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "HETATOM")):
            atom = parse_pdb_line(line)
            chain, orig_res = atom["chain_id"], atom["residue_number"]

            if chain not in chain_map: chain_map[chain] = {"last": orig_res, "new": start_index}
            elif chain_map[chain]["last"] != orig_res:
                chain_map[chain]["new"] += 1
                chain_map[chain]["last"] = orig_res

            atom["residue_number"] = chain_map[chain]["new"]

            new_lines.append(unparse_pdb_line({
                "record": atom["record"], "serial": atom["atom_number"], "name": atom["atom_name"],
                "altLoc": "", "resName": atom["residue_name"], "chainID": atom["chain_id"], "resSeq": atom["residue_number"],
                "x": atom["x"], "y": atom["y"], "z": atom["z"], "occupancy": atom["occupancy"], "tempFactor": atom["temp_factor"], "element": atom["element"]
            }))
        else: new_lines.append(line)
    return new_lines

# Obtener seq proteica (por cadenas)
def get_sequence(structure, ignore_chains=[]):
    seqs = {}
    for model in structure:
        for chain in model:
            if chain.id in ignore_chains: continue
            seqs[chain.id] = ""
            for res in chain:
                if not res.resname in residues_names: continue
                for i in range(res.id[1]-1-len(seqs[chain.id])): seqs[chain.id] += "X"  # Rellenar gaps con X                
                seqs[chain.id] += aminotools.convert_amino_acid_letter(res.resname)
    return seqs

# Escribir archivo .PDB a partir de una estructura
def write_structure_to_pdb(structure, outfile):
    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile)
    print(f"PDB file successfully generated at {outfile}")

# Mergear múltiples estructuras de Bio.PDB (structures, i.e. list) en una sola
def merge_structures(structures):
    builder = StructureBuilder.StructureBuilder()
    builder.init_structure("Merged")
    builder.init_model(0)
    
    chain_id = ord("A")
    for i, structure in enumerate(structures):
        for model in structure:
            for chain in model:
                new_id = chr(chain_id)
                builder.init_chain(new_id)
                for residue in chain:
                    builder.structure[0][new_id].add(residue.copy())
                chain_id += 1

    merged_structure = builder.get_structure()
    
    io = PDBIO()
    io.set_structure(merged_structure)

    return merged_structure

# Convertir .pdb a .cif
def pdb_to_cif(input_pdb, outdir=None):
    outfile = input_pdb.replace(".pdb", ".cif")
    if not outdir is None: outfile = os.path.join(outdir, os.path.basename(outfile))
    structure = get_structure(input_pdb)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(outfile)
    print(f"PDB succesfully converted to CIF at {outfile}")

# Convertir de .cif a .pdb
def cif_to_pdb(input_cif, outdir=None):
    outfile = input_cif.replace(".cif", ".pdb")
    if not outdir is None: outfile = os.path.join(outdir, os.path.basename(outfile))
    structure = get_structure(input_cif, format="cif")
    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile)
    print(f"CIF successfully converted to PDB at {outfile}")

####################################################################################### PARSING PDBS ##################################################################################################

# Obtener estructura de PDBParser
def get_structure(pdb_path, format="pdb"):
    if format == "pdb": parser = PDBParser(QUIET=True)  # QUIET=True evita warnings innecesarios
    elif format == "cif": parser = MMCIFParser(QUIET=True)
    else: raise ValueError(f"format '{format}' is not a valid format. Try whether with '.pdb' or '.cif' format")

    structure = parser.get_structure("model", pdb_path)
    return structure

# Mapea una lista de residuos (en formato Position.Chain e.g. '147C') sobre una estructura y devuelve los objetos de los residuos
def map(structure, residues_list):
    return [res for res in get_residues(structure) if f"{res.id[1]}{res.parent.id}" in residues_list]

# Medir distancia entre dos residuos
def distance(res1, res2):
    return np.linalg.norm(res1.center_of_mass() - res2.center_of_mass())

# Medir distancias
def calculate_distance(coords1, coords2):
    try: return np.linalg.norm(coords1 - coords2)
    except: raise ValueError(f"coords must be a 3D list. {coords1} vs {coords2}")

# Obtener residuos de una estructura PDB
def get_residues(structure):
    residues = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.resname in residues_names: residues.append(res)
    return residues

# Obtener heteroatomos de una estructura PDB
def get_non_residues(structure, ignore=["HOH"]):
    non_residues = []
    for model in structure:
        for chain in model:
            for res in chain:
                if not res.resname in residues_names: 
                    if not res.resname in ignore: non_residues.append(res)
    return non_residues

# Obtiene los C-a
def get_CA(structure):
    CA = []
    residues = get_residues(structure)
    for res in residues:
        for atom in res:
            if atom.id == "CA": CA.append(atom)
    return CA

# Obtener los atomos con un dado ID de los residuos de la proteina (si ids esta vacia incluye todos los atomos). as_dict lo devuelve en forma de diccionario, sino es una lista
def get_atoms(structure, ids=[], include_hetatoms=False, as_dict=False):
    atoms = {} if as_dict else []
    
    # Obtener residuos estándar
    residues = get_residues(structure)
    
    # Agregar HETATM si se solicita
    if include_hetatoms: residues += get_non_residues(structure)
    
    for res in residues:
        res_id = res.id[1]
        chain_id = res.parent.id
        full_id = f"{chain_id}{res_id}"
        
        # Filtrado por ID si corresponde
        if ids and res_id not in ids: continue
        
        if as_dict:
            if full_id not in atoms: atoms[full_id] = []
        else: atom_list = []
        
        for atom in res:
            if as_dict: atoms[full_id].append(atom)
            else: atom_list.append(atom)
        
        if not as_dict and atom_list: atoms.append(atom_list)
    
    return atoms

# Obtener coordenadas de lista de atomos
def get_coords(atoms):
    return [atom.coord for atom in atoms]

# Calcular RMSD entre 2 listas de coordenadas
def calculate_rmsd(coords1, coords2):
    coords1 = np.asarray(coords1)
    coords2 = np.asarray(coords2)
    assert coords1.shape == coords2.shape, "Coordinate lists must have the same shape"
    diff = coords1 - coords2
    return np.sqrt((diff ** 2).sum() / len(coords1))

# RMSD para estructuras desorientadas
def kabsch_rmsd(P, Q):
    P = np.array(P); Q = np.array(Q)
    assert P.shape == Q.shape, "Coordinate arrays must match"

    # Center both point clouds
    P_centered = P - P.mean(axis=0); Q_centered = Q - Q.mean(axis=0)

    # Compute covariance matrix
    C = np.dot(P_centered.T, Q_centered)

    # Singular Value Decomposition
    V, S, Wt = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0
    if d: V[:, -1] *= -1

    # Rotation matrix U
    U = np.dot(V, Wt)

    # Rotate P
    P_rot = np.dot(P_centered, U)

    # Compute RMSD
    rmsd = np.sqrt(np.mean(np.sum((P_rot - Q_centered)**2, axis=1)))
    return rmsd

# Alineamiento estructural
def align_structures(mobile_structure, ref_atoms, mobile_atoms):
    if len(ref_atoms) != len(mobile_atoms):
        raise ValueError("Atom lists must have the same length.")

    # Use the Atom objects directly
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, mobile_atoms)

    # Deep copy to preserve original mobile structure
    aligned_mobile_structure = copy.deepcopy(mobile_structure)

    # Apply transformation to all atoms in the copied structure
    for model in aligned_mobile_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(super_imposer.rotran[0], super_imposer.rotran[1])

    return aligned_mobile_structure

# Obtener residuos cercanos a uno o mas ligando/s
def get_active_site(structure, ligands_names=[], threshold=8, include_distance=False, include_resname=False):
    ligands = [res for res in get_non_residues(structure) if res.resname in ligands_names]
    residues = get_residues(structure)
    active_site = {}
    for res in residues:
        for ligand in ligands:
            dist = distance(res, ligand)
            if dist <= threshold: 
                res_id = f"{res.resname}{res.id[1]}{res.parent.id}"
                lig_id = f"{ligand.resname}{ligand.id[1]}{ligand.parent.id}"
                if not include_resname: 
                    res_id = res_id.replace(res.resname, "")
                    lig_id = lig_id.replace(ligand.resname, "")
                if include_distance: res_id = [res_id, dist]
                
                if not lig_id in active_site: active_site[lig_id] = [res_id]
                else: active_site[lig_id].append(res_id)    
    return active_site

# Obtener vecinos de una lista de residuos 'target_residues' (la lista es en formato POSITION.CHAIN e.g. 170A)
def get_neighbors(structure, target_residues, threshold=8, include_distance=False, include_resname=False):
    residues = get_residues(structure)
    non_residues = get_non_residues(structure)
    _target_residues = [res for res_list in [residues, non_residues] for res in res_list if f"{res.id[1]}{res.parent.id}" in target_residues]
    neighbors = {}
    for target_res in _target_residues:
        for res_list in [residues, non_residues]:
            for res in res_list:
                if f"{res.id[1]}{res.parent.id}" == f"{target_res.id[1]}{target_res.parent.id}": continue
                dist = distance(target_res, res)
                if dist <= threshold:
                    key = f"{target_res.resname}{target_res.id[1]}{target_res.parent.id}"
                    val = f"{res.resname}{res.id[1]}{res.parent.id}"
                    if not include_resname: key = key.replace(target_res.resname, ""); val.replace(res.resname, "")
                    if include_distance: val = [val, dist]
                    if key not in neighbors: neighbors[key] = [val]
                    else: neighbors[key].append(val)
    return neighbors

# Obtener residuos de interfaz
def get_interface(structure, threshold=8, ignore_chains=[], include_distance=False, include_resname=False):
    residues = get_residues(structure)
    chains = []
    residues_by_chain = []
    for res in residues:
        if res.parent.id in ignore_chains: continue
        if res.parent.id not in chains:
            chains.append(res.parent.id)
            residues_by_chain.append([])
        residues_by_chain[-1].append(res)

    interface = {}
    for i, chain in enumerate(residues_by_chain):
        for res in chain:
            for j, chain2 in enumerate(residues_by_chain):
                if i == j: continue
                for res2 in chain2:
                    dist = distance(res, res2)
                    if dist <= threshold:
                        key = f"{res.resname}{res.id[1]}{res.parent.id}"
                        val = f"{res2.resname}{res2.id[1]}{res2.parent.id}"
                        if not include_resname: key = key.replace(res.resname, ""); val.replace(res2.resname, "")
                        if include_distance: val = [val, dist]
                        if key not in interface: interface[key] = [val]
                        else: interface[key].append(val)
    return interface

# Obtener residuos de la superficie
def get_surface(structure, threshold=5, include_distance=False, include_resname=False):
    surface_matrix = GetRDSurface(structure[0])  # Get Residue Depth's Surface Model
    residues = get_residues(structure)           # Get Structure's Residues
    surface = []
    for res in residues:
        dist = min_dist(res.center_of_mass(), surface_matrix)
        if dist <= threshold: 
            res_id = f"{res.resname}{res.id[1]}{res.parent.id}"
            if not include_resname:  res_id = res_id.replace(res.resname, "")
            if include_distance: res_id = [res_id, dist]
            surface.append(res_id)  
    return surface

# Obtener residuos segun su grupo quimico (chem_class). main_chemical es un bool que indica si queres filtrar segun su grupo quimico principal o no (ej: PHE es principalmente aromatic pero tambien es hydrophobic)
def get_residues_by_chemical(residues, chem_class="aromatic", main_chemical=False):
    if not main_chemical: return [res for res in residues if chem_class.lower() in aminotools.get_all_chemical_classifications(res.resname)]
    else: return [res for res in residues if aminotools.get_main_chemical_classification(res.resname) == chem_class.lower()]

###################### CONTACTOS ################################

# Atomos terminales de cada residuo
aminoacids_terminal_atoms_dict = {
    "ALA": {"CB"}, "VAL": {"CG1", "CG2"}, "LEU": {"CD1", "CD2"}, "ILE": {"CD1"}, "MET": {"CE"}, "PRO": {"CD"},
    "PHE": {"CZ"}, "TYR": {"CZ", "OH"}, "TRP": {"CH2"}, "CYS": {"SG"},
    "SER": {"OG"}, "THR": {"OG1", "CG2"}, "ASN": {"OD1", "ND2"}, "GLN": {"OE1", "NE2"},
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "HIS": {"NE2", "CE1"}, "LYS": {"NZ"}, "ARG": {"NH1", "NH2"}, "GLY": set()
}

# Lista con estado de valencia de atomos
valence_dict = {
    "H": 1, "C": 4, "N": 3, "O": 2, "F": 1, "P": 3, "S": 2, "CL": 1,
    "BR": 1, "I": 1, "NA": 1, "K": 1, "MG": 2, "CA": 2, "ZN": 2, "FE": 2, "CU": 2
}

# Funcion para obtener atomos terminales de una dada lista de atomos
def get_terminal_atoms(atoms, valence_dict, threshold=2):
    terminal_atoms = set()
    grouped_atoms = {}
    for atom in atoms:
        res_id = (atom.parent.resname, atom.parent.id[1], atom.parent.parent.id)
        grouped_atoms.setdefault(res_id, []).append(atom)

    for res_id, res_atoms in grouped_atoms.items():
        coords = [a.coord for a in res_atoms if a.name[0] != "H"]
        tree = cKDTree(coords)
        for i, atom in enumerate(res_atoms):
            if atom.name[0] == "H": continue
            try:
                max_n = valence_dict.get(atom.name[:2], valence_dict.get(atom.name[0], 4))
            except: continue
            neighbors = tree.query_ball_point(atom.coord, threshold)
            if len(neighbors) - 1 < max_n:
                terminal_atoms.add(atom)
    return terminal_atoms

# Obtener interacciones hidrofobicas (C-C, C-S o S-S) entre 2 grupos de atomos
def get_hydrophobic_contacts(
    atoms1, atoms2, threshold=4.5,
    include_resname=False, include_distance=False, include_atom=False,
    terminal_atom_dict=aminoacids_terminal_atoms_dict,
    valence_dict=valence_dict
):
    contacts = {}

    # Precomputar terminales
    term_atoms1 = set()
    term_atoms2 = set()

    for atom in atoms1:
        if atom.parent.resname in terminal_atom_dict:
            if atom.name in terminal_atom_dict[atom.parent.resname] and atom.name[0] in {"C", "S"}:
                term_atoms1.add(atom)

    for atom in atoms2:
        if atom.parent.resname in terminal_atom_dict:
            if atom.name in terminal_atom_dict[atom.parent.resname] and atom.name[0] in {"C", "S"}:
                term_atoms2.add(atom)

    # Agregar terminales dinámicos si no están en el diccionario
    extra1 = get_terminal_atoms(atoms1, valence_dict)
    extra2 = get_terminal_atoms(atoms2, valence_dict)

    term_atoms1 |= {a for a in extra1 if a.name[0] in {"C", "S"}}
    term_atoms2 |= {a for a in extra2 if a.name[0] in {"C", "S"}}

    for atom1 in term_atoms1:
        if atom1.parent.resname == "GLY": continue
        for atom2 in term_atoms2:
            if atom2.parent.resname == "GLY": continue
            dist = calculate_distance(atom1.coord, atom2.coord)
            if dist > threshold: continue

            key = f"{atom1.parent.resname}{atom1.parent.id[1]}{atom1.parent.parent.id}" if include_resname else f"{atom1.parent.id[1]}{atom1.parent.parent.id}"
            value = f"{atom2.parent.resname}{atom2.parent.id[1]}{atom2.parent.parent.id}" if include_resname else f"{atom2.parent.id[1]}{atom2.parent.parent.id}"

            if include_atom:
                key = f"{key}_{atom1.name}"
                value = f"{value}_{atom2.name}"

            if include_distance:
                value = [value, dist]

            contacts.setdefault(key, []).append(value)

    return contacts

# Átomos aromáticos conocidos por residuo
aromatic_atoms_dict = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"}
    # Para ligandos, se detectará dinámicamente
}

# Para ligandos, encontrar atomos aromaticos
def get_aromatic_atoms_by_residue(residues):
    # Detección naive: 6 átomos planos conectados en ciclo con distancias ~1.4 Å
    # Aquí simplemente marcamos todos los C/N planos con 2 o 3 conexiones como candidatos
    aromatic_atoms = set()
    for residue in residues:
        atoms = [a for a in residue.get_atoms() if a.name[0] not in {"H"}]
        tree = cKDTree([a.coord for a in atoms])
        for i, atom in enumerate(atoms):
            neighbors = tree.query_ball_point(atom.coord, 1.6)
            if len(neighbors) in {2, 3} and atom.name[0] in {"C", "N"}:
                aromatic_atoms.add(atom)
    return aromatic_atoms

# Encontrar contactos aromaticos entre distintas listas de atomos
def get_aromatic_contacts(
    atoms1, atoms2, threshold=5.5,
    include_resname=False, include_distance=False, include_atom=False,
    aromatic_atom_dict=aromatic_atoms_dict
):
    contacts = {}

    def is_aromatic(atom):
        resname = atom.parent.resname
        if resname in aromatic_atom_dict:
            return atom.name in aromatic_atom_dict[resname]
        return atom in aromatic_dynamic_atoms

    # Precalcular átomos aromáticos dinámicos
    residues1 = {a.parent for a in atoms1 if a.parent.resname not in aromatic_atom_dict}
    residues2 = {a.parent for a in atoms2 if a.parent.resname not in aromatic_atom_dict}
    aromatic_dynamic_atoms = get_aromatic_atoms_by_residue(residues1 | residues2)

    aromatic1 = [a for a in atoms1 if is_aromatic(a)]
    aromatic2 = [a for a in atoms2 if is_aromatic(a)]

    for atom1 in aromatic1:
        for atom2 in aromatic2:
            if atom1 == atom2: continue
            dist = calculate_distance(atom1.coord, atom2.coord)
            if dist > threshold: continue

            key = f"{atom1.parent.resname}{atom1.parent.id[1]}{atom1.parent.parent.id}" if include_resname else f"{atom1.parent.id[1]}{atom1.parent.parent.id}"
            value = f"{atom2.parent.resname}{atom2.parent.id[1]}{atom2.parent.parent.id}" if include_resname else f"{atom2.parent.id[1]}{atom2.parent.parent.id}"

            if include_atom:
                key = f"{key}_{atom1.name}"
                value = f"{value}_{atom2.name}"

            if include_distance:
                value = [value, dist]

            contacts.setdefault(key, []).append(value)

    return contacts

# Obtener los puentes de hidrogeno entre dos listas de atomos
def get_hbond_contacts(
    atoms1, atoms2, threshold=3.5, min_angle=120, max_angle=180, donors_elements=("N", "O", "F", "S"),
    include_resname=False, include_distance=False, include_angle=False, include_atom=False
):
    contacts = {}

    def is_donor(atom): return (atom.name[0] in donors_elements) and ("H" in atom.name.upper())

    def get_angle(donor, acceptor):
        neighbors = [a for a in donor.parent.get_atoms() if a != donor]
        if neighbors:
            v1 = acceptor.coord - donor.coord
            v2 = neighbors[0].coord - donor.coord
            cosine = np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2))
            return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
        else: return 180  # asumo ideal si no hay vecinos

    def process_pair(donor, acceptor):
        dist = calculate_distance(donor.coord, acceptor.coord)
        if dist > threshold:
            return

        angle = get_angle(donor, acceptor)
        if not (min_angle <= angle <= max_angle):
            return

        key = f"{donor.parent.resname}{donor.parent.id[1]}{donor.parent.parent.id}" if include_resname else f"{donor.parent.id[1]}{donor.parent.parent.id}"
        value = f"{acceptor.parent.resname}{acceptor.parent.id[1]}{acceptor.parent.parent.id}" if include_resname else f"{acceptor.parent.id[1]}{acceptor.parent.parent.id}"

        if include_atom:
            key = f"{key}_{donor.name}"
            value = f"{value}_{acceptor.name}"

        if include_distance or include_angle:
            extras = []
            if include_distance:
                extras.append(dist)
            if include_angle:
                extras.append(angle)
            value = [value] + extras

        contacts.setdefault(key, []).append(value)

    for atom1 in atoms1:
        for atom2 in atoms2:
            if atom1 == atom2:
                continue
            # Caso 1: atom1 donante, atom2 aceptor
            if is_donor(atom1) and (atom2.name[0] in donors_elements):
                process_pair(atom1, atom2)
            # Caso 2: atom2 donante, atom1 aceptor
            if is_donor(atom2) and (atom1.name[0] in donors_elements):
                process_pair(atom2, atom1)

    return contacts