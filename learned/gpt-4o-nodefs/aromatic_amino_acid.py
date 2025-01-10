"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is characterized by the presence of an amino group, 
    a carboxylic acid group, and an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule has features of an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns to detect amino and carboxylic acid groups
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0+0]")  # NH2 or NH group, including ionized
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H0,H1-]")  # COOH, COO- or ester-like group
    aromatic_pattern = Chem.MolFromSmarts("a")  # Any aromatic atom (generic)

    # Check for the presence of an amino group
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Check for the presence of a carboxylic acid group
    # This pattern now includes ester or related groups to address false negatives
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or related group found"

    # Check for the presence of an aromatic ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"

    return True, "Contains amino group, carboxylic acid (or related), and aromatic ring"

# Testing with a variety of example SMILES auto-recognition
examples = [
    "Nc1ccc(cn1)C(O)=O", # 6-aminonicotinic acid
    "Nc1cc(O)cc(c1)C(O)=O", # 3-amino-5-hydroxybenzoic acid
    "N[C@@H](Cc1c[nH]cn1)C(O)=O", # L-histidine
    "CCCCNc1ccc(cc1)C(O)=O", # 4-(butylamino)benzoic acid
    "Nc1ccccc1C(O)=O", # anthranilic acid
]

for example in examples:
    print(is_aromatic_amino_acid(example))