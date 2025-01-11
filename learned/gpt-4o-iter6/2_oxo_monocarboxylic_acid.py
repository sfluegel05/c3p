"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid must have a monocarboxylic acid group and a 2-oxo substituent,
    with flexibility in how they are connected.
    
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
    
    # Check for a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for 2-oxo group on a carbon adjacent to the carboxyl carbon
    # Pattern: carboxylic acid (C(=O)O) adjacent to a carbon with a double bonded oxygen (oxo group)
    oxo_adjacent_pattern = Chem.MolFromSmarts("C(=O)C(=O)")
    if mol.HasSubstructMatch(oxo_adjacent_pattern):
        return True, "Contains carboxylic acid group with adjacent 2-oxo substituent"
    else:
        return False, "No proximate 2-oxo group to carboxylic acid group found"
        
# Example Testing
example_smiles = [
    "OC(=O)C(=O)CC1=CNC2=CC=C(O)C=C12",  # 3-(5-hydroxyindol-3-yl)pyruvic acid
    "CCCCCCCCC(=O)C(O)=O",               # 2-oxooctanoic acid
    "CC(C)(O)C(=O)C(O)=O",               # 3-hydroxy-3-methyl-2-oxopentanoic acid
    "OC(=O)C(=O)CP(O)(O)=O"              # 3-phosphonopyruvic acid
]

for smi in example_smiles:
    result, reason = is_2_oxo_monocarboxylic_acid(smi)
    print(f"SMILES: {smi} -> is_2_oxo_monocarboxylic_acid? {result}. Reason: {reason}")