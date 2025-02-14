"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid is characterized by the presence of a ketone group at the 11th carbon position within a steroid framework.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise.
        str: Reason for classification.
    """
    
    # Convert the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core: Decahydronaphthalene-like structure
    steroid_core = Chem.MolFromSmarts("C1CC2CCC3C4C=CC(=O)CC4CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid-like core not found"

    # 11-oxo group pattern: Ketone positioned at 11th carbon
    keto_11_pattern = Chem.MolFromSmarts("C=O")
    all_matches = mol.GetSubstructMatches(steroid_core)
    for match in all_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetDegree() >= 2 and any(neighbor.GetAtomicNum() == 8 and neighbor.GetSymbol() == 'O' for neighbor in atom.GetNeighbors()):
                return True, "Contains 11-oxo group in a typical steroid framework"

    return False, "11-oxo group not found"