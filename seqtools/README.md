
# Seqtools - Herramientas para Manejo de Secuencias FASTA y Alineamientos

**Seqtools** es una colección de funciones en Python para procesar archivos de secuencia, especialmente en formato FASTA.  
Permite leer, escribir, combinar archivos FASTA, realizar alineamientos múltiples con MAFFT y extraer residuos alineados.

---

## FUNCIONALIDADES PRINCIPALES

- Leer secuencias desde archivos FASTA a un diccionario.
- Escribir diccionarios `{id: secuencia}` como archivos FASTA.
- Combinar múltiples archivos FASTA en uno.
- Ejecutar MAFFT para generar alineamientos múltiples de secuencias (MSA).
- Obtener residuos perfectamente alineados desde un archivo MSA.

---

## USO BÁSICO

### Importación

```python
from seqtools import *
```

---

## EJEMPLOS DE USO

### 1. Leer un archivo FASTA

```python
seqs = open_fasta("secuencias.fasta")
print(seqs["seq1"])  # Muestra la secuencia asociada al ID 'seq1'
```

### 2. Escribir un FASTA a partir de un diccionario

```python
seqs_dict = {"id1": "MKT...", "id2": "GAV..."}
write_fasta(seqs_dict, "out.fasta")
```

### 3. Combinar múltiples archivos FASTA en uno solo

```python
merge_fastas(["a.fasta", "b.fasta", "c.fasta"], "combinado.fasta")
```

### 4. Generar un alineamiento múltiple con MAFFT

```python
mafft_MSA("combinado.fasta", "alineado.fasta")
# Requiere que MAFFT esté disponible en el PATH del sistema
```

### 5. Extraer residuos alineados (idénticos en todas las secuencias del MSA)

```python
aligned_residues = get_aligned_residues("alineado.fasta", add_pos=True, ref_id=0)
print(aligned_residues)  # Ejemplo: ['A5', 'L6', 'G7']
```

---

## LISTA DE FUNCIONES

- `open_fasta(fasta_file) → dict`
- `write_fasta(seqs_dict, outfile)`
- `merge_fastas(fasta_files, outfile)`
- `mafft_MSA(fasta_file, outfile)`
- `get_aligned_residues(msa_file, add_pos=True, ref_id=0) → list`

---

## REQUISITOS

- Python 3.x  
- MAFFT debe estar instalado y disponible en la línea de comandos para `mafft_MSA`.

---

## AUTOR

Franco Salvatore

## LICENCIA

Libre uso con atribución
