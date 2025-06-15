# This script/library consists on helping with getting any physicochemical property of aminoacids and other tools related to aminoacids

import math

# Converts one letter aminoacid to three letters and viceversa
def convert_amino_acid_letter(code):
    one_to_three = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }
    three_to_one = {v.lower(): k for k, v in one_to_three.items()}
    code = code.strip().capitalize()
    
    if len(code) == 1: return one_to_three.get(code.upper(), 'unknown').upper()
    elif len(code) == 3: return three_to_one.get(code.lower(), 'unknown').upper()
    else: return 'unknown'

# Gets the height and width of a given amino acid (1-letter code)
def get_dimensions(aa):
    aa_dimensions =  {
        'A': {'height': 5.7, 'width': 3.4}, 'R': {'height': 10.0, 'width': 4.8}, 'N': {'height': 6.0, 'width': 3.9}, 'D': {'height': 5.6, 'width': 4.2},
        'C': {'height': 5.5, 'width': 3.5}, 'E': {'height': 6.3, 'width': 4.1}, 'Q': {'height': 6.7, 'width': 4.1}, 'G': {'height': 4.5, 'width': 3.4},
        'H': {'height': 6.6, 'width': 4.2}, 'I': {'height': 6.0, 'width': 4.8}, 'L': {'height': 6.3, 'width': 4.2}, 'K': {'height': 9.8, 'width': 4.6},
        'M': {'height': 6.6, 'width': 4.2}, 'F': {'height': 8.0, 'width': 4.5}, 'P': {'height': 5.5, 'width': 4.2}, 'S': {'height': 5.2, 'width': 3.6},
        'T': {'height': 5.4, 'width': 3.8}, 'W': {'height': 9.4, 'width': 6.0}, 'Y': {'height': 8.1, 'width': 4.8}, 'V': {'height': 5.6, 'width': 4.2},
    }
    return aa_dimensions.get(aa.upper(), None)

# Gets the volume of a given amino acid (1-letter code) assuming it as a sphere
def get_volume(aa):
    dimensions = get_dimensions(aa)
    if not dimensions or 'height' not in dimensions or 'width' not in dimensions: return None
    diameter = (dimensions['height'] + dimensions['width']) / 2
    radius = diameter / 2
    volume = (4 / 3) * math.pi * (radius ** 3)
    return volume

# Determines the chemical group(s) of a given amino acid (1-letter code)

def get_main_chemical_classification(aa):
    return {
        'D': 'negative', 'E': 'negative',
        'K': 'positive', 'R': 'positive', 'H': 'positive',
        'S': 'hydrophilic', 'T': 'hydrophilic', 'N': 'hydrophilic',
        'Q': 'hydrophilic', 'C': 'hydrophilic', 'Y': 'aromatic',
        'F': 'aromatic', 'W': 'aromatic',
        'A': 'hydrophobic', 'I': 'hydrophobic', 'L': 'hydrophobic',
        'M': 'hydrophobic', 'V': 'hydrophobic',
        'G': 'hydrophobic', 'P': 'hydrophobic'
    }.get(aa.upper(), 'unknown')


def get_all_chemical_classifications(aa):
    aa = aa.upper()
    if len(aa) == 3: aa = convert_amino_acid_letter(aa)
    groups = {
        'positive': {'K', 'R', 'H'},
        'negative': {'D', 'E'},
        'hydrophilic': {'C', 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Y', 'W'},
        'hydrophobic': {'A', 'C', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'P', 'V', 'W', 'Y'},
        'aromatic': {'F', 'H', 'W', 'Y'},
        'aliphatic': {'A', 'I', 'L', 'P', 'V'}
    }
    return [g for g, residues in groups.items() if aa in residues] or ['unknown']

# Returns the average mass of a given amino acid (1-letter code)
def get_mass(aa):
    avg_masses = {
        'A': 71.0788,  'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
        'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
        'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
        'M': 131.1926, 'F': 147.1766, 'P': 97.1167,  'S': 87.0782,
        'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }
    return avg_masses.get(aa.upper(), 'unknown')

# Returns the Kyte-Doolittle hydropathy index of a given amino acid (1-letter code)
def get_hydropathy_index(aa):
    hydropathy = {
        'I': 4.5,  'V': 4.2,  'L': 3.8,  'F': 2.8,  'C': 2.5,
        'M': 1.9,  'A': 1.8,  'G': -0.4, 'T': -0.7, 'S': -0.8,
        'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
        'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }
    return hydropathy.get(aa.upper(), 'unknown')

# Returns the corresponding codon for a given amino acid (1-letter code) or viceversa
def map_codon(query):
    # Diccionario: aminoácido -> codones (ADN)
    aa_to_codon = {
        'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'],
        '*': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']
        }

    # Invertir el diccionario: codón -> aminoácido
    codon_to_aa = {codon: aa for aa, codons in aa_to_codon.items() for codon in codons}

    query = query.upper().replace('U', 'T')  # Convertir ARN a ADN si es necesario

    if len(query) == 1:  # aminoácido
        return aa_to_codon.get(query, f"Amino acid '{query}' not found.")
    elif len(query) == 3:  # codón
        return codon_to_aa.get(query, f"Codon '{query}' not found.")
    else:
        return "Input must be either a single-letter amino acid or a 3-letter codon."

# Gives a list of possible mRNA mutations for a given input codon and alternative aminoacid of interest
def get_possible_mutations(codon, aa):
    codon = codon.upper()
    aa = aa.upper()
    
    if len(codon) != 3 or aa not in 'ACDEFGHIKLMNPQRSTVWY':
        return "Invalid input. Codon must be 3 letters and amino acid must be a valid one-letter code."
    
    # Get the list of codons for the target amino acid
    target_codons = map_codon(aa)
    
    mutations = []
    for i in range(3):  # Iterate over each position in the codon
        for base in 'ACGT':
            if base != codon[i]:  # Change only if it's different
                new_codon = codon[:i] + base + codon[i+1:]
                if new_codon in target_codons:
                    mutations.append(new_codon)
    
    return mutations