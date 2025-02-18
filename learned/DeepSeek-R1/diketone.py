"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains exactly two ketone functionalities.
    A ketone is defined as a carbonyl group (C=O) where the carbonyl carbon is bonded to two other carbon atoms.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a ketone group: carbon bonded to two other carbons and double-bonded to oxygen
    ketone_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])=O")
    matches = mol.GetSubstructMatches(ketone_pattern)
    ketone_count = len(matches)
    
    if ketone_count == 2:
        return True, "Contains exactly two ketone groups"
    else:
        return False, f"Found {ketone_count} ketone groups, expected 2"