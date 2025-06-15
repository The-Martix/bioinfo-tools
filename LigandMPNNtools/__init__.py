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

## ANALISIS ##

# Esta funcion filtra las pred_seqs de una dada run que no cumplan con un valor minimo de una dada metric
def filter_pred_seqs(pred_seqs, metric="ovr_conf", threshold=0.5):
    filtered_pred_seqs = []
    for pred_seq in pred_seqs:
        value = getattr(pred_seq, metric)
        if value < threshold: continue
        else: filtered_pred_seqs.append(pred_seq)
    return filtered_pred_seqs

# Esta funcion hace un scatterplot de las metricas elegidas para una dada run
def plot_run(run, metrics=["ovr_conf", "seq_conf"]):
    x_values = [getattr(pred_seq, metrics[0]) for pred_seq in run.pred_seqs]
    y_values = [getattr(pred_seq, metrics[1]) for pred_seq in run.pred_seqs]
    labels = [f"run_{run.id}_{pred_seq.id}" for pred_seq in run.pred_seqs]
    plt.scatter(x_values, y_values)

    # Agregar etiquetas a cada punto
    for i, label in enumerate(labels): plt.text(x_values[i], y_values[i], label, fontsize=8)
    plt.xlabel(metrics[0]); plt.ylabel(metrics[1]); plt.title(f"Run {run.id}")
    plt.tight_layout();  plt.grid(True)
    plt.show()

# Esta funcion hace un scatterplot de las metricas elegidas para una lista de runs (los thresholds se usan para filtrar por valores de las metricas correspondientes)
def plot_runs(runs, metrics=["ovr_conf", "seq_conf"], thresholds=[0, 0]):
    x_vals, y_vals, labels = [], [], []
    for run in runs:
        pred_seqs = run.pred_seqs
        for metric, threshold in zip(metrics, thresholds): pred_seqs = filter_pred_seqs(pred_seqs, metric, threshold)
        x_values = [getattr(pred_seq, metrics[0]) for pred_seq in pred_seqs]
        y_values = [getattr(pred_seq, metrics[1]) for pred_seq in pred_seqs]
        labs = [f"run_{run.id}_{pred_seq.id}" for pred_seq in pred_seqs]
        x_vals.extend(x_values); y_vals.extend(y_values); labels.extend(labs)

    plt.scatter(x_vals, y_vals)
    # Agregar etiquetas a cada punto
    for i, label in enumerate(labels): plt.text(x_vals[i], y_vals[i], label, fontsize=8)
    plt.xlabel(metrics[0]); plt.ylabel(metrics[1])
    plt.tight_layout();  plt.grid(True)
    plt.show()

# Esta funcion grafica la confianza por posicion de residuo de una dada seq predicha por LigandMPNN
def plot_conf_per_res(pred_seq):
    # Extraer datos
    confs = [res.conf for res in pred_seq.residues]
    poses = [res.pos for res in pred_seq.residues]

    # Estilo
    sns.set(style="whitegrid")

    # Crear figura
    plt.figure(figsize=(14, 6))
    scatter = plt.scatter(range(len(poses)), confs, c=confs, cmap='viridis', s=80, edgecolors='k', alpha=0.9)

    # Agregar número de residuo encima de cada punto
    for i, (x, y) in enumerate(zip(range(len(poses)), confs)):  plt.text(x, y + 0.02, str(poses[i]), ha='center', va='bottom', fontsize=8)

    # Estética de ejes
    plt.xlabel("Residue Index");  plt.ylabel("Confidence per Position"); plt.title("Predicted Confidence by Residue")
    plt.colorbar(scatter, label="Confidence Score")
    plt.xticks([]); plt.ylim(-0.1, 1.1); plt.tight_layout()
    plt.show()

# Esta funcion grafica la confianza promedio por posicion de todas las seqs predichas en una run de LigandMPNN
def plot_mean_conf_per_res(seqs):
    # Extraer datos
    poses = [res.pos for res in seqs[0].residues]
    confs = []
    for i in range(len(seqs[0].residues)):
        _confs = []
        for seq in seqs: _confs.append(seq.residues[i].conf)
        confs.append(mean(_confs))
        
    # Estilo
    sns.set(style="whitegrid")

    # Crear figura
    plt.figure(figsize=(14, 6))
    scatter = plt.scatter(range(len(poses)), confs, c=confs, cmap='viridis', s=80, edgecolors='k', alpha=0.9)

    # Agregar número de residuo encima de cada punto
    for i, (x, y) in enumerate(zip(range(len(poses)), confs)):  plt.text(x, y + 0.02, str(poses[i]), ha='center', va='bottom', fontsize=8)

    # Estética de ejes
    plt.xlabel("Residue Index"); plt.ylabel("Mean Confidence per Position"); plt.title(f"Predicted Mean Confidence by Residue")
    plt.colorbar(scatter, label="Confidence Score")
    plt.xticks([]); plt.ylim(-0.1, 1.1);  plt.tight_layout()
    plt.show()

# Esta funcion obtiene un diccionario de los residuos por posicion y las confianzas medias para cada aminoacido disenado por LigandMPNN en cada posicion 
def get_residues_per_pos(seqs, positions=[]):
    len_seq = len(seqs[0].seq)
    if len(positions) == 0: positions = list(range(1, len_seq+1))  # Si no se especifican posiciones, usar todas las posiciones de la secuencia
    residues_per_pos = {}

    for pos in range(1, len_seq+1):
        if not pos in positions: continue
        pos_dict = {}

        for seq in seqs:
            residue = seq.residues[pos-1]
            res_name = residue.name
            pos_dict.setdefault(res_name, []).append(residue.conf)

        # Promediar y ordenar por confianza descendente
        averaged = {res_name: mean(confs) for res_name, confs in pos_dict.items()}
        sorted_avg = dict(sorted(averaged.items(), key=lambda item: item[1], reverse=True))
        residues_per_pos[pos] = sorted_avg

    return residues_per_pos

def get_redesigned_residues_proportions(seqs, positions=[]):
    len_seq = len(seqs[0].seq)
    if len(positions) == 0: positions = list(range(1, len_seq+1))  # Si no se especifican posiciones, usar todas las posiciones de la secuencia
    residues_per_pos = {}

    def get_residue_proportion(residues):
        proportions = {}
        for res in set(residues):
            proportion = residues.count(res) / len(residues)
            proportions[res] = proportion
        return proportions

    for pos in range(1, len_seq+1):
        if not pos in positions: continue
        residues = []
        for seq in seqs: residues.append(seq.residues[pos-1].name)
        residues_per_pos[pos] = get_residue_proportion(residues)

    return residues_per_pos

# Esta funcion grafica un barplot que representa las confianzas por posicion promedio de una run de LigandMPNN
def barplot_conf_per_res(seqs, seq_range="ALL", cmap="YlOrRd"):
    residues_per_pos = get_residues_per_pos(seqs, seq_range)
    sns.set(style="whitegrid")
    n_positions = len(residues_per_pos)
    fig, ax = plt.subplots(figsize=(max(10, n_positions * 0.8), 6))

    # Preparar datos
    positions, aas, confs = [], [], []
    for pos, res_dict in residues_per_pos.items():
        for aa, conf in res_dict.items(): positions.append(pos); aas.append(aa); confs.append(conf)

    # Escalado a [0, 1] para el colormap
    norm = plt.Normalize(0, 1)
    colormap = cm.get_cmap(cmap)

    # Coordenadas x para cada posición
    unique_positions = sorted(residues_per_pos.keys())
    pos_indices = {pos: i for i, pos in enumerate(unique_positions)}
    x = [pos_indices[p] for p in positions]

    # Asignar color basado en confianza
    bar_colors = [colormap(norm(conf)) for conf in confs]
    bars = ax.bar(x, confs, color=bar_colors, edgecolor='black', alpha=0.95)
    # Etiquetas arriba de cada barra
    for bar, aa in zip(bars, aas):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 0.02, aa, ha='center', va='bottom', fontsize=8)

    # Ejes y título
    ax.set_xticks(list(pos_indices.values()))
    ax.set_xticklabels([str(pos) for pos in unique_positions])
    ax.set_xlabel("Residue Position"); ax.set_ylabel("Confidence"); ax.set_title("Confidence per Residue and Mutation (Warmmap)")
    ax.set_ylim(0, 1.05)

    # Barra de colores (opcional)
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm);  sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax); cbar.set_label("Confidence")

    plt.tight_layout(); plt.show()
