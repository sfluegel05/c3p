"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid must have a carboxylic acid group and a 2-oxo substituent
    at the adjacent position.

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

    # Define SMILES patterns for carboxylic acid and 2-oxo adjacent to carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    oxo_adjacent_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)")
    
    # Check for carboxylic acid
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for 2-oxo adjacent to carboxylic acid
    if not mol.HasSubstructMatch(oxo_adjacent_pattern):
        return False, "No 2-oxo group adjacent to carboxylic group found"
    
    return True, "Contains carboxylic acid group and 2-oxo group at adjacent position"

# Examples
example_smiles = [
    "OC(=O)C(=O)CC1=CNC2=CC=C(O)C=C12",  # 3-(5-hydroxyindol-3-yl)pyruvic acid
    "CCCCCCCCC(=O)C(O)=O"               # 2-oxooctanoic acid
]

for smi in example_smiles:
    result, reason = is_2_oxo_monocarboxylic_acid(smi)
    print(f"SMILES: {smi} -> is_2_oxo_monocarboxylic_acid? {result}. Reason: {reason}")