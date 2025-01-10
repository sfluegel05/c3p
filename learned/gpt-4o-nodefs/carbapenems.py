"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem typically has a β-lactam fused to a five-membered ring and often a sulfur atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for carbapenem: a β-lactam fused to a five-membered ring
    # Here, the SMARTS pattern looks for a 4-membered ring around a nitrogen (the β-lactam)
    # fused to another ring, which is the common feature in carbapenems
    carbapenem_pattern = Chem.MolFromSmarts("C1C2N1C(=O)[CX3]2")  # Simplified pattern (could be refined)
    
    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "Does not match the carbapenem structural motif"

    # Check for sulfur presence if necessary
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if sulfur_count == 0:
        return False, "No sulfur atom found, less likely to be a typical carbapenem"

    return True, "Contains the β-lactam core with a fused five-membered ring typical of carbapenems"

# Note: The SMARTS pattern and logic can be further refined to capture more specific features of carbapenems