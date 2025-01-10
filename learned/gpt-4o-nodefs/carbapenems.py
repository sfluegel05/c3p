"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem typically features a bicyclic ring structure with a β-lactam
    fused to a five-membered ring.

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
    
    # Define a refined SMARTS pattern for carbapenem
    # This pattern identifies a β-lactam (4-membered) fused with a five-membered ring
    carbapenem_pattern = Chem.MolFromSmarts("C1C(=O)N2C[C@H](C)C12")  # Broad example pattern
    
    if not mol.HasSubstructMatch(carbapenem_pattern):
        return False, "Does not match the carbapenem structural core motif"

    # Check for optional sulfur
    has_sulfur = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
    if not has_sulfur:
        return True, "Matches carbapenem core structure, but lacks typical sulfur"

    return True, "Matches carbapenem core structure with potential sulfur"

# Further refinements can be made based on a more specific understanding
# and evaluation against known carbapenem structures.