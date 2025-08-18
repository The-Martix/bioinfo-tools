# This script parse AlphaFold3 (AF3) runs and its metrics. By Franco Salvatore, 2025

import os
import structools
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import defaultdict
import math

def load_json(file):
    f = open(file, "r")
    data = json.load(f)
    f.close()
    return data

class Run:
    def __init__(self, pdb_path, _id, model):
        self.pdb_path = pdb_path
        self.id = _id
        self.model = model
        self.chains = []
        self.confidences_path = self.pdb_path.replace(f"model_{model}.pdb", f"summary_confidences_{model}.json")
        self.get_confidences()

    def get_confidences(self):
        self.conf_data = load_json(self.confidences_path)
        self.iptm = self.conf_data['iptm']
        self.pTM = self.conf_data['ptm']

    def parse_structure(self, correct_resname={"B": {"LG1": "HEM"}, "C": {"LG1": "EST"}}):
        """
        Parseo usando structools.open_pdb y structools.parse_pdb_line.
        Se evita Bio.PDB para no colapsar HETATM con nombres repetidos.
        """
        lines = structools.open_pdb(self.pdb_path)

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
            plddt = d["temp_factor"]  # B-factor como proxy de plddt o score
            coords = [d["x"], d["y"], d["z"]]

            residue.add_atom(unique_atom_name, atom_serial, plddt, coords)

        # calcular promedios al final
        for ch in self.chains:
            for r in ch.residues:
                r.get_mean_plddt()
            ch.get_mean_plddt()
        self.mean_plddt = float(np.mean([ch.mean_plddt for ch in self.chains]))

class Chain:
    def __init__(self, _id):
        self.id = _id
        self.residues = []
    
    def get_mean_plddt(self):
        if self.residues:
            self.mean_plddt = float(np.mean([res.mean_plddt for res in self.residues]))
        else:
            self.mean_plddt = float("nan")

class Residue:
    def __init__(self, resname, respos):
        self.resname = resname
        self.respos = respos
        self.atoms = []
        self.mean_plddt = None

    def add_atom(self, atomname, atomid, plddt, coords):
        self.atoms.append(Atom(atomname, atomid, plddt, coords))

    def get_mean_plddt(self):
        if self.atoms:
            self.mean_plddt = float(np.mean([atom.plddt for atom in self.atoms]))
        else:
            self.mean_plddt = float("nan")

class Atom:
    def __init__(self, atomname, atomid, plddt, coords):
        self.name = atomname
        self.id = atomid
        self.plddt = plddt
        self.coords = coords

# Parse run
def parse_run(folder_path, correct_resname={"B": {"LG1": "HEM"}, "C": {"LG1": "EST"}}, model=0, keep_pdb=True):
    '''
    Parsea una corrida de AF3. Generalmente la corrida devuelve una carpeta (folder_path) con los modelos que genera y sus confianzas y metricas
    '''

    cif_file = [f for f in os.listdir(folder_path) if f"model_{model}.cif" in f][0]
    cif_path = os.path.join(folder_path, cif_file)

    # Generate pdb files from cif files
    pdb_path = cif_path.replace(".cif", ".pdb")
    structools.cif_to_pdb(cif_path, show_print=False)

    run = Run(pdb_path, os.path.splitext(os.path.basename(pdb_path))[0], model)
    run.parse_structure(correct_resname=correct_resname)

    if not keep_pdb:
        os.remove(pdb_path)
    return run

# Plot plDDT
def plot_plddt(run,
              level="atom",                 # "atom" | "residue" | "chain"
              color_by_chain=True,          # sólo aplica cuando level="atom"
              show_chain_mean=True,          # sólo aplica cuando level="atom"
              mean_color="black",            # color de la línea de la media por cadena
              cmap_name="tab10",             # colormap para cadenas (cuando color_by_chain=True)
              marker_size=12,                # tamaño de marcador
              line_alpha=0.9,
              fig=None, ax=None,
              show=True):
    """
    Plotea plddt de un objeto Run.
    - level="atom": scatter de átomos (opcional: colorear por cadena y trazar medias por cadena con labels).
    - level="residue": un punto por residuo (mean plddt del residuo).
    - level="chain": barras por cadena (mean plddt de la cadena).
    """

    if run is None or not hasattr(run, "chains") or len(run.chains) == 0:
        raise ValueError("The 'run' has no chains to plot.")

    # Crear fig/ax si hace falta
    if ax is None or fig is None:
        fig, ax = plt.subplots(figsize=(10, 4 if level != "chain" else 5))

    # Recolectar datos con orden estable (como quedaron en run.chains)
    # También guardamos offsets para saber los rangos de cada cadena en el eje X (útil para labels)
    x_vals = []
    y_vals = []
    chain_ids = []
    chain_ranges = {}  # chain_id -> (x_start, x_end) inclusive
    x_cursor = 0

    if level == "atom":
        for ch in run.chains:
            start = x_cursor
            for r in ch.residues:
                for a in r.atoms:
                    if a.plddt is None or (isinstance(a.plddt, float) and math.isnan(a.plddt)):
                        continue
                    x_vals.append(x_cursor)
                    y_vals.append(float(a.plddt))
                    chain_ids.append(ch.id)
                    x_cursor += 1
            end = x_cursor - 1
            if end >= start:
                chain_ranges[ch.id] = (start, end)

        if len(x_vals) == 0:
            raise ValueError("There are no plddt values per atom.")

        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)
        chain_ids = np.array(chain_ids)

        if color_by_chain:
            # Mapeo de cadena -> color usando colormap
            unique_chains = list(chain_ranges.keys())
            cmap = cm.get_cmap(cmap_name, max(1, len(unique_chains)))
            color_map = {cid: cmap(i) for i, cid in enumerate(unique_chains)}
            colors = [color_map[cid] for cid in chain_ids]
            ax.scatter(x_vals, y_vals, s=marker_size, c=colors, alpha=line_alpha, edgecolors="none")

            # Medias por cadena (a nivel atómico)
            if show_chain_mean:
                for i, cid in enumerate(unique_chains):
                    rng = chain_ranges[cid]
                    mask = (x_vals >= rng[0]) & (x_vals <= rng[1])
                    if np.any(mask):
                        mean_val = float(np.mean(y_vals[mask]))
                        # línea horizontal de la media
                        ax.hlines(mean_val, xmin=rng[0], xmax=rng[1],
                                  colors=mean_color, linestyles="--", linewidth=1.5)
                        # label cerca del extremo derecho del rango
                        ax.text(rng[1] + max(2, int(0.01 * len(x_vals))), mean_val,
                                f"{cid} mean={mean_val:.2f}", va="center", fontsize=9, color=mean_color)
            # Marcar separadores entre cadenas
            for cid, rng in chain_ranges.items():
                ax.axvline(rng[0] - 0.5, color="lightgray", linewidth=0.8, linestyle=":")
            # Leyenda manual con puntos de ejemplo por cadena
            # (Opcional; muchos puntos pueden saturar la leyenda)
            # Construimos handles únicos por cadena
            handles = []
            labels = []
            for cid in unique_chains:
                handles.append(plt.Line2D([0], [0], marker='o', linestyle='',
                                          markersize=marker_size/3, color=color_map[cid]))
                labels.append(f"Chain {cid}")
            ax.legend(handles, labels, title="Chains", frameon=False, fontsize=9)
        else:
            # Sin color por cadena: un solo color por default
            ax.scatter(x_vals, y_vals, s=marker_size, alpha=line_alpha, edgecolors="none")

        ax.set_xlabel("Atoms (ordered by atomid)")
        ax.set_ylabel("plddt (B-factor)")
        ax.set_title(f"plddt at atomic level — run {run.id}")

        # Ajustes de límites y estética
        ax.set_xlim(-1, max(x_vals) + 10)
        # Intento de rango Y razonable si todo es igual
        if np.allclose(np.nanmin(y_vals), np.nanmax(y_vals)):
            ax.set_ylim(np.nanmin(y_vals) - 1, np.nanmax(y_vals) + 1)

    elif level == "residue":
        # Un punto por residuo con su mean_plddt
        for ch in run.chains:
            start = x_cursor
            for r in ch.residues:
                if r.mean_plddt is None or (isinstance(r.mean_plddt, float) and math.isnan(r.mean_plddt)):
                    continue
                x_vals.append(x_cursor)
                y_vals.append(float(r.mean_plddt))
                chain_ids.append(ch.id)
                x_cursor += 1
            end = x_cursor - 1
            if end >= start:
                chain_ranges[ch.id] = (start, end)

        if len(x_vals) == 0:
            raise ValueError("There are no plddt values per residue.")

        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)

        ax.plot(x_vals, y_vals, marker="o", linestyle="-", alpha=line_alpha)
        for cid, rng in chain_ranges.items():
            ax.axvline(rng[0] - 0.5, color="lightgray", linewidth=0.8, linestyle=":")

        ax.set_xlabel("Residues (ordered by resid)")
        ax.set_ylabel("mean plddt by residue")
        ax.set_title(f"plddt at residual level — run {run.id}")

        if np.allclose(np.nanmin(y_vals), np.nanmax(y_vals)):
            ax.set_ylim(np.nanmin(y_vals) - 1, np.nanmax(y_vals) + 1)

    elif level == "chain":
        # Una barra por cadena con su mean_plddt (ya calculada en run.parse_structure)
        labels = []
        means = []
        for ch in run.chains:
            labels.append(str(ch.id))
            # Si querés asegurar media a nivel atómico en lugar de la media de residuos:
            # atomic_vals = [a.plddt for r in ch.residues for a in r.atoms if a.plddt is not None and not math.isnan(a.plddt)]
            # means.append(float(np.mean(atomic_vals)) if len(atomic_vals) else np.nan)
            means.append(float(ch.mean_plddt) if ch.residues else np.nan)

        y_vals = np.array(means, dtype=float)
        x = np.arange(len(labels))
        ax.bar(x, y_vals, alpha=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlabel("Chains")
        ax.set_ylabel("mean plddt by chain")
        ax.set_title(f"plddt at chain level — run {run.id}")
        for xi, yi in zip(x, y_vals):
            if not (isinstance(yi, float) and math.isnan(yi)):
                ax.text(xi, yi, f"{yi:.2f}", ha="center", va="bottom", fontsize=9)

        # Rango Y si todo igual
        if np.allclose(np.nanmin(y_vals), np.nanmax(y_vals)):
            ax.set_ylim(np.nanmin(y_vals) - 1, np.nanmax(y_vals) + 1)

    else:
        raise ValueError("level must be 'atom', 'residue' o 'chain'.")

    ax.set_ylim(bottom=0, top=105)
    ax.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax