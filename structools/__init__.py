# This script/library gets the PDB structure of a given pdb file and set structural classifications measuring distances of CAMs within a given threshold

from Bio.PDB import PDBParser
from Bio.PDB import Superimposer
from Bio.PDB import PDBIO
from Bio.PDB.ResidueDepth import get_surface as GetRDSurface
from Bio.PDB.ResidueDepth import min_dist
import copy
import aminotools
import numpy as np

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
def write_pdb(pdb_file, new_lines):
    file = open(pdb_file, 'w')
    for line in new_lines: file.write(line)
    file.close()
    print(f"Se ha escrito el archivo {pdb_file} con exito.")

# Obtengo las lineas de atomos del ligando del Docking
def get_ligand_atoms(lines):
    atom_lines = []
    for line in lines:
        if line.startswith("ATOM"):
            new_line = line.replace("ATOM", "HETATM")
            atom_lines.append(new_line)
    return atom_lines

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
def merge_ligand_to_pdb(reference_pdb, ligand_pdb):
    atom_lines = get_ligand_atoms(open_pdb(ligand_pdb))
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
def get_sequence(structure):
    seqs = {}
    for model in structure:
        for chain in model:
            seqs[chain.id] = ""
            for res in chain:
                if not res.resname in residues_names: continue
                seqs[chain.id] += aminotools.convert_amino_acid_letter(res.resname)
    return seqs

# Escribir archivo .PDB a partir de una estructura
def write_structure_to_pdb(structure, outfile):
    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile)
    print(f"PDB file successfully generated at {outfile}")

####################################################################################### PARSING PDBS ##################################################################################################

# Obtener estructura de PDBParser
def get_structure(pdb_path):
    parser = PDBParser(QUIET=True)  # QUIET=True evita warnings innecesarios
    structure = parser.get_structure("model", pdb_path)
    return structure

# Mapea una lista de residuos (en formato Position.Chain e.g. '147C') sobre una estructura y devuelve los objetos de los residuos
def map(structure, residues_list):
    return [res for res in get_residues(structure) if f"{res.id[1]}{res.parent.id}" in residues_list]

# Medir distancia entre dos residuos
def distance(res1, res2):
    return np.linalg.norm(res1.center_of_mass() - res2.center_of_mass())

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

# Obtener los atomos con un dado ID de los residuos de la proteina (si ids esta vacia incluye todos los atomos)
def get_atoms(structure, ids=[], include_hetatoms=False):
    atoms = {}
    residues = get_residues(structure)
    if include_hetatoms:
        non_residues = get_non_residues(structure)
        for res in non_residues: residues.append(res)
    for res in residues:
        for atom in res:
            if atom.id in ids or len(ids) == 0:
                try: atoms[f"{res.parent.id}{res.id[1]}"].append(atom)
                except: atoms[f"{res.parent.id}{res.id[1]}"] = [atom]
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

# Obtener vecinos de una lista de residuos (la lista es en formato POSITION.CHAIN e.g. 170A)
def get_neighbors(structure, target_residues, threshold=8, include_distance=False, include_resname=False):
    residues = get_residues(structure)
    _target_residues = [res for res in residues if f"{res.id[1]}{res.parent.id}" in target_residues]
    neighbors = {}
    for target_res in _target_residues:
        for res in residues:
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