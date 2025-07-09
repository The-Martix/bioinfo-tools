# ChemTools - Utilidades para Manejo de Compuestos Químicos y Archivos PDB/SDF

Este paquete reúne varias herramientas útiles para trabajar con archivos químicos en formatos PDB y SDF.  
Incluye funciones para extraer ligandos específicos de archivos PDB, convertir ligandos a SMILES, y combinar múltiples archivos SDF en uno solo.
En el futuro se espera seguir contribuyendo con esta libreria para hacerla mas completa. Por ejemplo, que se pueda convertir smiles a sdf y viceversa, etc 

---

## FUNCIONALIDADES PRINCIPALES

- Extraer ligandos específicos desde archivos `.pdb`.
- Guardar solo los ligandos seleccionados en un nuevo archivo `.pdb`.
- Convertir un archivo `.pdb` de un ligando a su representación SMILES.
- Combinar múltiples archivos `.sdf` en un único archivo de salida.
- Asegurar la correcta separación de moléculas en archivos SDF con la marca `$$$$`.

---

## USO BÁSICO

### Importación de funciones

```python
from chemtools import write_ligands_pdb, pdb_to_smiles, merge_sdf
```

### 1. Extraer ligandos de interés desde un archivo PDB

```python
write_ligands_pdb("complejo.pdb", ["HEM", "ATP"], "ligandos.pdb")
```

### 2. Convertir un ligando en formato PDB a SMILES

```python
smiles = pdb_to_smiles("ligando.pdb")
print(smiles)
```

### 3. Combinar múltiples archivos SDF

```python
merge_sdf(["lig1.sdf", "lig2.sdf"], "combinado.sdf")
```

---

## DEFINICIÓN DE FUNCIONES

### `write_ligands_pdb(input_pdb, ligands: list[str], output_pdb)`

Extrae solo las líneas correspondientes a los ligandos deseados desde un archivo `.pdb` y las guarda en otro archivo `.pdb`.  
Utiliza funciones de `structools`.

### `pdb_to_smiles(input_pdb: str, removeHs=False) -> str`

Convierte un archivo `.pdb` de un ligando a su representación SMILES utilizando RDKit.  
Lanza un error si el archivo no se puede cargar correctamente.

### `merge_sdf(input_files: list[str], output_file: str)`

Concatena múltiples archivos `.sdf` en un único archivo de salida.  
Asegura que cada molécula esté separada por `$$$$`.

```python
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
```

---

## REQUISITOS

- Python 3.9 o superior  
- [RDKit](https://www.rdkit.org/)  
- `structools`: módulo propio o externo para manejo de archivos PDB

---

## AUTOR

Franco Salvatore

## LICENCIA

Libre uso con atribución
