# Structools - Herramientas para Parsing, Análisis y Edición de Archivos PDB

**Structools** es una colección de funciones para facilitar el análisis estructural de proteínas en formato PDB.  
Permite extraer residuos, medir distancias, alinear estructuras, identificar sitios activos e interfaces,  
así como editar directamente archivos PDB para modificar, filtrar o reescribir estructuras.

---

## 🧰 Requisitos

- Python >= 3.x  
- Biopython  
- NumPy  
- ResidueDepth de Biopython (MSMS)  
- [`aminotools`](https://github.com/SalvaFran/bioinfo-tools/blob/main/structools/README.md) (ya sea en el `PATH` o en la misma carpeta)

**Instalación de dependencias:**

```bash
pip install biopython numpy
```

---

## 🔍 ¿Qué se puede hacer con Structools?

- Leer estructuras PDB.
- Extraer residuos, heteroátomos, átomos específicos o coordenadas.
- Calcular distancias entre residuos.
- Calcular RMSD entre estructuras (con y sin alineamiento previo).
- Alinear estructuras completas.
- Obtener residuos cercanos a un ligando o de la interfaz entre cadenas.
- Detectar residuos en la superficie (requiere `ResidueDepth`).
- Editar archivos PDB: remover aguas, combinar estructuras, escribir nuevos PDB.

---

## 🧪 Ejemplos de uso

### Importación

```python
from structools import *

structure = get_structure("1abc.pdb")
```

### Parsing: Obtener residuos del sitio activo

```python
# Obtener residuos a 6 Å del ligando HEM
site = get_active_site(structure, ligands_names=["HEM"], threshold=6)
print(site)
```

### Parsing: Calcular RMSD entre dos estructuras

```python
structure1 = get_structure("wt.pdb")
structure2 = get_structure("mut.pdb")

ca1 = get_CA(structure1)
ca2 = get_CA(structure2)

coords1 = get_coords(ca1)
coords2 = get_coords(ca2)

rmsd = calculate_rmsd(coords1, coords2)
print("RMSD sin alinear:", rmsd)
```

### Parsing: Alinear dos estructuras y calcular RMSD

```python
aligned = align_structures(structure2, ca1, ca2)

ca2_aligned = get_CA(aligned)
coords2_aligned = get_coords(ca2_aligned)

rmsd_aligned = calculate_rmsd(coords1, coords2_aligned)
print("RMSD alineado:", rmsd_aligned)
```

### Edición: Remover aguas y guardar un nuevo PDB

```python
lines = read_pdb("1abc.pdb")
clean_lines = remove_water(lines)
write_pdb("1abc_no_water.pdb", clean_lines)
```

### Edición: Combinar dos estructuras en un mismo PDB

```python
lines1 = read_pdb("receptor.pdb")
lines2 = read_pdb("ligando.pdb")

merged_lines = merge_pdbs(lines1, lines2)
write_pdb("complejo.pdb", merged_lines)
```

---

## 🧩 Funciones disponibles

```python
get_structure(path)
get_residues(structure)
get_non_residues(structure)
get_CA(structure)
get_atoms(structure, ids=[], include_hetatoms=False)
get_coords(atoms)
calculate_rmsd(coords1, coords2)
kabsch_rmsd(P, Q)
align_structures(mobile_structure, ref_atoms, mobile_atoms)
get_active_site(structure, ligands_names, threshold)
get_neighbors(structure, target_residues, threshold)
get_interface(structure, threshold)
get_surface(structure)
read_pdb(path)
write_pdb(path, lines)
remove_water(lines)
merge_pdbs(lines1, lines2)
```

---

## 📄 Autor y licencia

**Autor:** Franco Salvatore  
**Licencia:** Libre uso con atribución
"""
