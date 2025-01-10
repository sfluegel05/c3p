"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is characterized by an oxo (ketone) group at the third position on the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid backbone structure
    # Steroid backbone: three six-membered rings followed by a five-membered ring
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4=CC(=O)CCC34")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # Check for the 3-oxo group
    # This checks for a carbonyl group (=O) bound to the third carbon which is a part of the steroid backbone
    oxo_group_pattern = Chem.MolFromSmarts("C2=COCC(=O)CC3C2CCC4")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found"

    return True, "Molecule is a 3-oxo steroid"

# Example usage: Check if a given SMILES string is a 3-oxo steroid
# Call the function with a SMILES string as argument, for example:
# is_3_oxo_steroid('CC(=O)[C@H]1CC[C@H]2[C@@H]3CC=C4C(F)(F)C(=O)CC[C@]4(C)[C@H]3CC[C@]12C')
# This would return (True, "Molecule is a 3-oxo steroid") if the SMILES string matches the criteria.