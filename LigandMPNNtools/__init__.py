# This script/library processes LigandMPNN runs

import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import aminotools

# Remueve elementos vacios de una lista
def remove_gaps(lista):
    new_list = []
    for elem in lista: 
        if elem != "" and elem != None: new_list.append(elem.replace("\n", "").strip())
    return new_list

def mean(lista): return sum(lista)/len(lista)

# Class

class Run:
    def __init__(self, id, dir):
        self.id = id
        self.dir = dir
        self.proteins = []

class Protein:
    def __init__(self, id):
        self.id = id
        self.seqs = []
        self.template_seq = []

    def add_seq(self, seq): self.seqs.append(seq)

    def sort_pred_seqs(self, by="ovr_conf"): return sorted(self.seqs, key=lambda seq: getattr(seq, by), reverse=True)

    def get_mean_conf(self, by="ovr_conf"): return mean([getattr(seq, by) for seq in self.seqs])

    # Esta funcion obtiene la seq con mayor confianza de una run, a partir de las seqs de esa proteina
    def generate_more_conf_seq(self):
        residues_per_pos = get_residues_per_pos(self.seqs)
        best_res_per_pos = {}
        for pos in list(residues_per_pos.keys()): best_res_per_pos[pos] = { list(residues_per_pos[pos].keys())[0]: list(residues_per_pos[pos].values())[0] }
        more_conf_seq = ""
        for res_dict in list(best_res_per_pos.values()): more_conf_seq += aminotools.convert_amino_acid_letter(list(res_dict.keys())[0])
        
        return best_res_per_pos, more_conf_seq   

class Seq:
    def __init__(self, id, name, dir, seq, ovr_conf, lig_conf, seq_rec):
        self.id = id
        self.name = name
        self.dir = dir
        self.seq = seq
        self.ovr_conf = ovr_conf
        self.lig_conf = lig_conf
        self.seq_rec = seq_rec
        self.residues = []
        self.atoms_conf, self.hetatoms_conf = None, None

    def calculate_identity(self, seq):
        # Calculate identity between self.seq and seq
        matches = sum(1 for a, b in zip(self.seq, seq) if a == b)
        return matches / len(seq) if len(seq) > 0 else 0
    
class Residue:
    def __init__(self, name, pos, chain, conf):
        self.name = name
        self.pos = pos
        self.chain = chain
        self.conf = conf


## PBD ##

# Esta funcion devuelve las lineas de un pdb
def open_design_pdb(pdb_file):
    with open(pdb_file, 'r') as file: lines = file.readlines()
    return lines

# Esta funcion toma las lineas de la seq proteica de un pdb predicho por LigandMPNN y calcula su confianza global
def get_pred_conf(lines):
    atom_confidences, hetatom_confidences = [], []
    residues = []
    last_res_confs = []
    last_res = ""
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            splits = remove_gaps(line.split(" "))
            try: conf = float(splits[10])
            except: conf = float(splits[9])

            if line.startswith("ATOM"):
                atom_confidences.append(conf)

                resname, reschain, respos = splits[3], splits[4], int(splits[5])
                res = f"{resname} {reschain} {respos}"
                
                if res != last_res: 
                    if len(last_res_confs) > 0: 
                        last_res_conf = mean(last_res_confs)
                        last_res_confs = []
                        resname, reschain, respos = last_res.split(" ")
                        residue = Residue(resname.strip(), int(respos.strip()), reschain.strip(), last_res_conf)
                        residues.append(residue)
                    last_res = res
                
                last_res_confs.append(conf)
            
            elif line.startswith("HETATM"): 
                
                if last_res != None:
                    last_res_conf = mean(last_res_confs)
                    resname, reschain, respos = last_res.split(" ")
                    residue = Residue(resname.strip(), int(respos.strip()), reschain.strip(), last_res_conf)
                    residues.append(residue)
                    last_res = None

                hetatom_confidences.append(conf)
    
    atom_conf = mean(atom_confidences)
    hetatom_conf = mean(hetatom_confidences)

    return atom_conf, hetatom_conf, residues 

## FASTA ##

# Esta funcion devuelve en un diccionario las seqs de un fasta con su header:seq
def parse_fasta(fasta_file):
    fasta_dict = {}
    with open(fasta_file, 'r') as file: lines = file.readlines()
    for line in lines:
        if line.startswith('>'):
            header = line.strip()
            fasta_dict[header] = ''
        else:
            fasta_dict[header] += line.strip()
    return fasta_dict

# Esra funcion parsea una secuencia fasta predicha y obtiene su metadata
def parse_seq(seq_key):
    splits = seq_key.split(',')
    _id = int(splits[1].split("=")[1].strip())
    name = splits[0][1:].strip()
    for i, split in enumerate(splits):
        if "overall_confidence" in split:
            ovr_conf = float(splits[i].split("=")[1].strip())
            lig_conf = float(splits[i+1].split("=")[1].strip())
            seq_rec = float(splits[i+2].split("=")[1].strip())
            break
    return _id, name, ovr_conf, lig_conf, seq_rec

# Esta funcion parsea la corrida entera de LigandMPNN
def parse_run(run_dir):
    run_id = os.path.basename(run_dir)           # Obtener el id de la corrida
    run = Run(id=run_id, dir=run_dir)            # Crear el objeto Run
    fasta_path = rf"{run_dir}\seqs"              # Obtener el path de la carpeta de fasta
    fasta_files = os.listdir(fasta_path)         # Obtener los nombres de los archivos fasta 
    for fasta_file in fasta_files:
        #run.structure_id = fasta_file.split(".")[0]  # Obtener el id de la estructura usada en el LigandMPNN
        protein = Protein(fasta_file.split(".")[0])
        run.proteins.append(protein)
        fasta_dict = parse_fasta(f"{fasta_path}/{fasta_file}")  # Parsear el fasta y obtener el diccionario con las secuencias y sus encabezados        
        # Agregar la secuencia input y las predichas
        template_seq = list(fasta_dict.values())[0]
        protein.template_seq = template_seq
        for i, key in enumerate(list(fasta_dict.keys())):
            if i == 0: continue
            seq = fasta_dict[key]
            _id, name, ovr_conf, lig_conf, seq_rec = parse_seq(key)
            pred_seq = Seq(id=_id, name=name, dir=f"{fasta_path}/{fasta_file}", seq=seq, ovr_conf=ovr_conf, lig_conf=lig_conf, seq_rec=seq_rec)
            pdb_file = rf"{run.dir}/backbones/{protein.id}_{pred_seq.id}.pdb"
            pred_seq.atoms_conf, pred_seq.hetatoms_conf, pred_seq.residues = get_pred_conf(open_design_pdb(pdb_file=pdb_file))
            protein.add_seq(pred_seq)

    return run