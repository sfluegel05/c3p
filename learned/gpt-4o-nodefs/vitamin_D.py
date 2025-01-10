"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for secosteroid pattern: secosteroid B-ring opening and triene system
    # A simple SMARTS pattern for vitamin D secosteroid backbone could be a broken B-ring (C-C-C in sterols),
    # and conjugated system (C=C-C=C-C) in the secosteroid backbone
    
    secosteroid_pattern = Chem.MolFromSmarts("C1=CC2CCCCC2C3(C=CC=CC3)C1")
    
    if mol.HasSubstructMatch(secosteroid_pattern):
        return True, "Contains secosteroid structural motif common to vitamin D"
    
    return False, "Does not contain secosteroid structural motif"

# Example usage
smiles_str = "[C@@H]1(C[C@@H](C/C(=C/C=C/2\CCC[C@]3([C@]2(CC[C@]3([H])[C@](CCCC4(O)CCCC4)([H])C)[H])C)/C1=C)O)O"
result, reason = is_vitamin_D(smiles_str)
print(f"Result: {result}, Reason: {reason}")