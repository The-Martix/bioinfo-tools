
# Aminotools - Utilidades para Propiedades de Aminoácidos y Mutaciones

**Aminotools** es una librería de funciones en Python diseñada para facilitar el acceso a propiedades fisicoquímicas de aminoácidos, conversiones entre formatos y herramientas básicas para análisis de codones y posibles mutaciones.

---

## FUNCIONALIDADES PRINCIPALES

- Conversión entre códigos de una y tres letras para aminoácidos.
- Obtener dimensiones (alto y ancho) y calcular volumen aproximado.
- Determinar clasificación química principal y todas las clases posibles.
- Obtener la masa promedio de un aminoácido.
- Calcular el índice de hidropatía (Kyte-Doolittle).
- Mapear aminoácidos a codones y viceversa.
- Obtener posibles mutaciones de mRNA que produzcan un cambio deseado.

---

## USO BÁSICO

### Importación

```python
from aminotools import *
```

---

## EJEMPLOS DE USO

### 1. Conversión de código

```python
convert_amino_acid_letter("A")      # → 'ALA'
convert_amino_acid_letter("tyr")    # → 'Y'
```

### 2. Dimensiones y volumen

```python
get_dimensions("W")                 # → {'height': 9.4, 'width': 6.0}
get_volume("W")                     # → volumen estimado en Å³
```

### 3. Clasificación química

```python
get_main_chemical_classification("K")   # → 'positive'
get_all_chemical_classifications("W")   # → ['hydrophilic', 'hydrophobic', 'aromatic']
```

### 4. Masa y propiedades hidrofóbicas

```python
get_mass("V")                      # → 99.1326
get_hydropathy_index("V")         # → 4.2
```

### 5. Mapeo de codones

```python
map_codon("R")                    # → ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
map_codon("TGG")                  # → 'W'
```

### 6. Mutaciones posibles

```python
get_possible_mutations("AAG", "E")
# → Lista de codones que con una mutación convierten AAG (Lys) en Glu (E)
```

---

## LISTA DE FUNCIONES

- `convert_amino_acid_letter(code)`
- `get_dimensions(aa)`
- `get_volume(aa)`
- `get_main_chemical_classification(aa)`
- `get_all_chemical_classifications(aa)`
- `get_mass(aa)`
- `get_hydropathy_index(aa)`
- `map_codon(query)`
- `get_possible_mutations(codon, aa)`

---

## AUTOR

Franco Salvatore

## LICENCIA

Libre uso con atribución
