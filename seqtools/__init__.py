# Este script es para realizar procesos relacionados con secuencias (e.g. fasta) como por ej alineamientos, etc

import os

# Extraer seq de un fasta
def open_fasta(fasta_file):
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
    os.system(f"mafft {fasta_file} > {outfile}")
    print(f"MSA succesfully generated at {outfile}")

# Obtener residuos que alinean (i.e. matchean) en un MSA. add_pos agrega la posicion del residuo alineado. red_id es la seq que usa como referencia del MSA
def get_aligned_residues(msa_file, add_pos=True, ref_id=0):
    seqs_dict = open_fasta(msa_file)
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

