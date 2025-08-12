# Script para parsear corridas de RoseTTAFold-All-Atom.
# Lee archivos PDB generados por el modelo, extrae información de cadenas, residuos y átomos,
# normaliza nombres de residuos según un mapeo definido, asigna identificadores únicos a átomos repetidos
# y calcula el lDDT medio por residuo a partir de los valores de B-factor.

import os
import structools
import numpy as np

class Run:
    def __init__(self, _id):
        self.id = _id
        self.chains = []
        self.pdb_path = None

    def parse_structure(self, pdb_path, correct_resname={"B": {"LG1": "HEM"}, "C": {"LG1": "EST"}}):
        """
        Parseo usando structools.open_pdb y structools.parse_pdb_line.
        Se evita Bio.PDB para no colapsar HETATM con nombres repetidos.
        """
        self.pdb_path = pdb_path
        lines = structools.open_pdb(pdb_path)

        # Mapa rápido para no duplicar objetos
        chains_map = {}  # chain_id -> Chain
        # (chain_id, resSeq, resNameNormalizado) -> Residue
        residues_map = {}
        # contadores por residuo para renombrar atomos repetidos (C1, C2...), por elemento
        elem_counters = {}  # key_res -> dict(element -> count)

        for line in lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            d = structools.parse_pdb_line(line)
            chain_id = d["chain_id"]
            resname_raw = d["residue_name"]
            resseq = d["residue_number"]

            # Normalización opcional de nombre de residuo por cadena
            resname = resname_raw
            if chain_id in correct_resname and resname_raw in correct_resname[chain_id]:
                resname = correct_resname[chain_id][resname_raw]

            # Crear/obtener Chain
            if chain_id not in chains_map:
                chains_map[chain_id] = Chain(chain_id)
                self.chains.append(chains_map[chain_id])

            # Crear/obtener Residue (clave incluye resSeq y resName normalizado)
            res_key = (chain_id, resseq, resname)
            if res_key not in residues_map:
                residue = Residue(resname=resname, respos=resseq)
                chains_map[chain_id].residues.append(residue)
                residues_map[res_key] = residue
                elem_counters[res_key] = {}

            residue = residues_map[res_key]

            # Elemento: del campo element; si falta, infiere de atom_name
            element = (d["element"] or d["atom_name"].strip()[:2]).strip()
            if len(element) == 0:
                element = "X"
            # normalizar a una o dos letras tipo PDB (primera mayúscula, segunda minúscula si existe)
            element = element[0].upper() + (element[1:].lower() if len(element) > 1 else "")

            # Contador por elemento dentro del residuo para crear nombres únicos
            cnt = elem_counters[res_key].get(element, 0) + 1
            elem_counters[res_key][element] = cnt
            unique_atom_name = f"{element}{cnt}"

            atom_serial = d["atom_number"]
            lddt = d["temp_factor"]  # B-factor como proxy de lDDT o score

            residue.add_atom(unique_atom_name, atom_serial, lddt)

        # calcular promedios al final
        for ch in self.chains:
            for r in ch.residues:
                r.get_mean_lddt()
            ch.get_mean_lddt()

class Chain:
    def __init__(self, _id):
        self.id = _id
        self.residues = []
    
    def get_mean_lddt(self):
        if self.residues:
            self.mean_lddt = float(np.mean([res.mean_lddt for res in self.residues]))
        else:
            self.mean_lddt = float("nan")

class Residue:
    def __init__(self, resname, respos):
        self.resname = resname
        self.respos = respos
        self.atoms = []
        self.mean_lddt = None

    def add_atom(self, atomname, atomid, lddt):
        self.atoms.append(Atom(atomname, atomid, lddt))

    def get_mean_lddt(self):
        if self.atoms:
            self.mean_lddt = float(np.mean([atom.lddt for atom in self.atoms]))
        else:
            self.mean_lddt = float("nan")

class Atom:
    def __init__(self, atomname, atomid, lddt):
        self.name = atomname
        self.id = atomid
        self.lddt = lddt

def parse_run(pdb_path, correct_resname={"B": {"LG1": "HEM"}, "C": {"LG1": "EST"}}):
    run = Run(os.path.splitext(os.path.basename(pdb_path))[0])
    run.parse_structure(pdb_path, correct_resname=correct_resname)
    return run

# Ejemplo de uso
# run = parse_run("path/a/pdb")