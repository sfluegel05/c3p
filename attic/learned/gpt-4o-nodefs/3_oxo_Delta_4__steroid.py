"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CC2=CC(C=C3C2CCC4C3(CCC4)[H])[C@](C1)([C@H]4[C@@](C)C(=O))")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for the 3-oxo group (C=O at position 3)
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C@]1([H])")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "Missing 3-oxo group"
    
    # Look for the Delta(4) double bond (C=C between C4 and C5)
    delta4_pattern = Chem.MolFromSmarts("C=CC")
    delta4_matches = mol.GetSubstructMatches(delta4_pattern)
    if not delta4_matches:
        return False, "Missing Delta(4) double bond"

    return True, "Contains the structural features of a 3-oxo-Delta(4) steroid"

# Example of how the program might be used:
# result, reason = is_3_oxo_Delta_4__steroid("[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])C(=O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO")
# print(result, reason)