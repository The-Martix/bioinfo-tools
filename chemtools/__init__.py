# Este script es para trabajar con archivos/compuesto quimicos

# imports
import structools
from rdkit import Chem
from rdkit.Chem import AllChem

# Esta funcion toma un PDB en el que estan tus ligandos de interes, una lista con los nombres de los ligandos y devuelve un pdb con solo los ligandos de interes
def write_ligands_pdb(input_pdb, ligands: list[str], output_pdb):
    lines = structools.open_pdb(input_pdb)
    ligand_lines = structools.get_hetatoms_lines(lines, hetatoms=ligands)
    structools.write_pdb(output_pdb, ligand_lines)

# Esta funcion toma un archivo PDB con SOLO los ligandos de interes y genera su SMILES correspondiente
def pdb_to_smiles(input_pdb: str, removeHs=False) -> str:
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=removeHs)
    if mol: return Chem.MolToSmiles(mol)
    else: raise ValueError("Error al cargar el PDB.")

# Esta funcion toma una lista con las rutas de los archivos SDF de entrada y los concatena en un solo archivo SDF de salida
def merge_sdf(input_files: list[str], output_file: str):
    with open(output_file, 'w') as fout:
        for file in input_files:
            with open(file, 'r') as fin:
                lines = fin.readlines()
                for line in lines:
                    fout.write(line)
                # Asegura que termina con salto de línea si no lo tiene
                if not lines[-1].strip() == "$$$$":
                    fout.write("$$$$\n")

    print(f"Archivos SDF combinados en {output_file}")

import subprocess

# Conversor de archivos moleculares
def convert_molecule_format(input_file, output_file):
    """
    Convierte un archivo molecular entre formatos PDB y SDF usando Open Babel CLI.
    
    Parámetros:
    - input_file (str): ruta al archivo de entrada (PDB o SDF).
    - output_file (str): ruta al archivo de salida (PDB o SDF).
    
    Requiere que 'obabel' esté instalado y en el PATH del sistema.
    """
    # Detectar formatos por extensión
    ext_in = input_file.split('.')[-1].lower()
    ext_out = output_file.split('.')[-1].lower()
    
    formats = {'pdb': 'pdb', 'sdf': 'sdf', 'mol2': 'mol2', 'mol': 'mol', 'smi': 'smi'}
    
    if ext_in not in formats or ext_out not in formats:
        raise ValueError("Formatos soportados: pdb, sdf, mol2, mol, smi")
    
    cmd = ['obabel', '-i' + formats[ext_in], input_file,
                  '-o' + formats[ext_out], '-O', output_file]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Archivo convertido: {input_file} → {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error en la conversión: {e}")
