"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid must have a carboxylic acid group and an oxo substituent
    adjacent to the carbonyl carbon of the acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMILES patterns for the functional groups
    # Carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    # 2-oxo pattern explicitly next to a carboxylic acid group
    oxo_adjacent_pattern = Chem.MolFromSmarts("C(=O)C(=O)")

    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for 2-oxo group adjacent to carboxylic acid
    if not mol.HasSubstructMatch(oxo_adjacent_pattern):
        return False, "No 2-oxo group adjacent to carboxylic group found"
    
    return True, "Contains carboxylic acid group and 2-oxo group adjacent to each other"

# Example Testing
example_smiles = [
    "OC(=O)C(=O)CC1=CNC2=CC=C(O)C=C12",  # 3-(5-hydroxyindol-3-yl)pyruvic acid
    "CCCCCCCCC(=O)C(O)=O",               # 2-oxooctanoic acid
    "CC(C)(O)C(=O)C(O)=O",               # 3-hydroxy-3-methyl-2-oxopentanoic acid
    "NC(=O)C(=O)C(O)=O",                 # oxaluric acid (incorrect, lacks oxo-adjacency pattern)
    "CCCCC(=O)C(O)=O",                   # 2-oxohexanoic acid
    "OC(=O)C(=O)CP(O)(O)=O"              # 3-phosphonopyruvic acid
]

for smi in example_smiles:
    result, reason = is_2_oxo_monocarboxylic_acid(smi)
    print(f"SMILES: {smi} -> is_2_oxo_monocarboxylic_acid? {result}. Reason: {reason}")