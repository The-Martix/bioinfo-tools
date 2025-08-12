# Este script es para realizar procesos relacionados con secuencias (e.g. fasta) como por ej alineamientos, etc

import os

# Lista de aminoacidos
amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X"]

# List de bases
dna_nucleotides = ["A", "T", "C", "U"]
rna_nucleotides = ["A", "U", "C", "G"]

# Esta funcion chequea que un archivo fasta sea valido (hay que indicar los tipos de seq del fasta en la variable 'seq_type' i.ie "protein", "dna", "rna")
def check_fasta(fasta_file, seq_type="protein", allow_gaps=False):
    """
    Verifica si un archivo FASTA contiene secuencias de aminoácidos o nucleótidos válidas.
    
    Parámetros:
    - fasta_file (str): Ruta al archivo FASTA.
    - seq_type (str): Tipo de secuencia esperada ("protein", "dna", "rna").
    
    Retorna:
    - bool: True si el archivo es válido, False en caso contrario.
    """
    valid_nucleotides = dna_nucleotides if seq_type == "dna" else rna_nucleotides if seq_type == "rna" else amino_acids
    with open(fasta_file, 'r') as file:
        for i, line in enumerate(file):
            if not line.startswith(">"):  # Ignorar líneas de encabezado
                for char in line.strip():
                    if char not in valid_nucleotides:
                        if char == "-" and allow_gaps: continue
                        return False, f"Invalid character '{char}' found in sequence '{line.strip()}' at line {i}. Expected characters: {valid_nucleotides}"
    return True, ""

# Esta funcion corrige errores comunes en archivos fasta, como espacios en blanco o caracteres no válidos.
def fix_fasta(fasta_file, seq_type="protein", outfile=None):
    if outfile is None: outfile = fasta_file
    """
    Corrige errores comunes en archivos FASTA, como espacios en blanco o caracteres no válidos.
    
    Parámetros:
    - fasta_file (str): Ruta al archivo FASTA.
    - seq_type (str): Tipo de secuencia esperada ("protein", "dna", "rna").
    
    Retorna:
    - str: Ruta al archivo FASTA corregido.
    """
    valid_chars = dna_nucleotides if seq_type == "dna" else rna_nucleotides if seq_type == "rna" else amino_acids
    new_lines = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                new_lines.append(line)
            else:
                cleaned = "".join([c.upper() if c.upper() in valid_chars else "X" for c in line])
                new_lines.append(cleaned)

    with open(outfile, 'w') as f:
        for line in new_lines:
            f.write(line + "\n")

    print(f"FASTA file fixed at {outfile}")

# Extraer seq de un fasta
def open_fasta(fasta_file, seq_type="protein", allow_gaps=False):
    check, msg = check_fasta(fasta_file, seq_type, allow_gaps=allow_gaps)
    if not check: raise ValueError(f"Invalid Fasta file. {msg}")
    seqs = {}
    file = open(fasta_file, "r")
    lines = file.readlines()
    file.close()
    last_id = ""
    for line in lines:
        if line.startswith(">"):
            last_id = line[1:].strip()
            seqs[last_id] = ""
        else:
            try: seqs[last_id] += line.strip()
            except: continue
    return seqs

# Escribir fasta a partir de diccionario con {id: seq}
def write_fasta(seqs_dict, outfile):
    file = open(outfile, "w")
    for _id, seq in zip(list(seqs_dict.keys()), list(seqs_dict.values())):
        file.write(f">{_id}\n")
        file.write(f"{seq}\n")
    file.close()
    print(f"fasta file succesfully generated at {outfile}")

# Juntar varios fastas (como listas) en uno
def merge_fastas(fasta_files, outfile):
    merged_file = open(outfile, "w")
    for file in fasta_files:
        fasta = open_fasta(file)
        for _id, seq in zip(list(fasta.keys()), list(fasta.values())):
            merged_file.write(f">{_id}\n")
            merged_file.write(f"{seq}\n")
    merged_file.close()
    print(f"Merged fasta succesfully generated at {outfile}")

# Genera un archivo fatsa con el MSA generado por mafft (mafft debe estar en el PATH)
def mafft_MSA(fasta_file, outfile):
    check, msg = check_fasta(fasta_file)
    if not check: raise ValueError(f"Invalid FASTA file. {msg}")
    os.system(f"mafft {fasta_file} > {outfile}")
    print(f"MSA succesfully generated at {outfile}")

# Obtener residuos que alinean (i.e. matchean) en un MSA. add_pos agrega la posicion del residuo alineado. red_id es la seq que usa como referencia del MSA
def get_aligned_residues(msa_file, add_pos=True, ref_id=0):
    seqs_dict = open_fasta(msa_file, allow_gaps=True)
    aligned_residues = []
    seqs = list(seqs_dict.values())
    for i in range(len(seqs[ref_id])):
        res_set = []
        for seq in seqs: res_set.append(seq[i])
        if len(list(set(res_set))) > 1: continue
        aligned_residues.append(seqs[ref_id][i])
        if add_pos: 
            pos = i + 1 - seqs[ref_id][:i].count("-")
            aligned_residues[-1] += f"{pos}"
    return aligned_residues

# Encuentra el residuo alineado para un un dado residuo de referencia
def find_aligned_residue(ref_seq, alt_seq, residue_number):
    """
    Encuentra con qué aminoácido/gap de la secuencia alternativa se alinea un residuo dado de la referencia.
    
    Parameters:
    - ref_seq (str): secuencia de referencia con gaps
    - alt_seq (str): secuencia alternativa con gaps (alineada a la referencia)
    - residue_number (int): número de residuo en la referencia (sin contar gaps)
    
    Returns:
    - tuple: (residuo_referencia, residuo_alternativa, posicion_alineada)
    """
    assert len(ref_seq) == len(alt_seq), "Las secuencias deben tener la misma longitud"
    
    count = 0  # contador de residuos sin gap en referencia
    for i, res in enumerate(ref_seq):
        if res != "-":
            count += 1
            if count == residue_number:
                return res, alt_seq[i], i+1  # residuo ref, residuo alt, posición en alineamiento
    
    return None  # si no se encuentra

