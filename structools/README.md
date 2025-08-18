# structools

`structools` es una biblioteca Python para el análisis, manipulación y clasificación de estructuras proteicas. Utiliza y extiende Biopython para ofrecer herramientas avanzadas en el tratamiento de archivos `.pdb` y `.cif`, integrando una arquitectura basada en programación orientada a objetos (OOP).

---

## Características Principales

- **Generador de Estructuras OOP**
  - Lectura y parsing de archivos `.pdb` en una jerarquía clara: `Structure > Chain > Residue > Atom`.
  - Cada nivel posee propiedades bioquímicas (masa, centro de masa, etc.) y permite manipulaciones locales y globales.

- **Edición y Procesamiento de Estructuras**
  - Fusión y división de cadenas.
  - Renombrado de ligandos y remoción de componentes (agua, HETATM, hidrógenos, etc).
  - Escritura a formatos `.pdb` o `.cif`.

- **Módulos de Clasificación Funcional**
  - `active_site`: identifica residuos cercanos a ligandos.
  - `interface`: detecta contactos inter-cadenas.
  - `surface`: clasifica residuos en la superficie mediante "residue depth".

**Nota:** Para estas funciones (`active_site`, `interface`, `surface`), la estructura de entrada debe estar generada mediante `PDBParser`.

- **Herramientas Analíticas**
  - Cálculo de distancias, RMSD, alineamientos estructurales.
  - Contactos hidrofóbicos, aromáticos y puentes de hidrógeno.

---

## Estructura OOP

La arquitectura OOP permite una representación jerárquica precisa:

### `Structure`
- Contiene múltiples `Chain`.
- Métodos:
  - `add_chain(chain)`: añade una cadena.
  - `merge_chains(chain_ids, new_chain_id, renumber)`: fusiona cadenas.
  - `split_chain(chain_id, residue_id_groups, ...)`: divide una cadena.
  - `rename_hetatoms(dict)`, `remove_waters()`, `remove_hydrogens()`, `remove_hetatoms()`, `remove_peptides()`: limpieza estructural.
  - `write_pdb(path)`: exporta archivo `.pdb`.

### `Chain`
- Contiene una lista de `Residue`.
- Propiedades:
  - `mass`, `center_of_mass`.

### `Residue`
- Contiene una lista de `Atom`.
- Propiedades:
  - `mass`, `center_of_mass`.

### `Atom`
- Posee propiedades como coordenadas, tipo de elemento, ocupancia, etc.
- Propiedad: `mass`.

Ejemplo de creación:
```python
from structools import get_structure
structure = get_structure("example.pdb", src="OOP")
print(structure.mass)
```

---

## Métodos Generales

- `get_structure(path, format="pdb", src="OOP"|"PDBParser")`: retorna una estructura.
- `generate_structure(data, id)`: genera estructura OOP desde diccionario.
- `pdb_to_cif(path)`, `cif_to_pdb(path)`: conversiones entre formatos.
- `write_structure_to_pdb(structure, path)`: guarda estructura.
- `align_structures(mobile, ref_atoms, mobile_atoms)`: alinea dos estructuras.
- `calculate_rmsd(coords1, coords2)`, `kabsch_rmsd(P, Q)`: cálculo de RMSD.
- `get_residues`, `get_non_residues`, `get_atoms`, `get_CA`: funciones auxiliares de extracción.

---

## Ejemplos de Uso

### Identificación de Sitio Activo (via `PDBParser`)
```python
from structools import get_structure, get_active_site
structure = get_structure("example.pdb", src="PDBParser")
active = get_active_site(structure, ligands_names=["ATP"])
```

### Residuo de Interfaz (via `PDBParser`)
```python
from structools import get_interface
interface = get_interface(structure)
```

### Residuos en Superficie (via `PDBParser`)
```python
from structools import get_surface
surface = get_surface(structure)
```

---

## Instalación
```bash
pip install biopython
# Agregar structools al PYTHONPATH o empaquetar como módulo
```

---

## Dependencias
- [Biopython](https://biopython.org/)
- NumPy
- SciPy

---

## Licencia
MIT

---

## Autor
Desarrollado por [Franco Salvatore] para aplicaciones de bioinformática estructural.

---

## Notas Finales
Este paquete está diseñado para usuarios avanzados que requieren herramientas especializadas en análisis estructural. Es especialmente útil para clasificación de ligandos, detección de sitios funcionales y modelado molecular.
