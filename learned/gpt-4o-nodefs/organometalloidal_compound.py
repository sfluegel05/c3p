"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound involves a metalloid element, we focus on arsenic (As) for this class.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of Arsenic (As)
    has_arsenic = any(atom.GetAtomicNum() == 33 for atom in mol.GetAtoms())
    if not has_arsenic:
        return False, "No arsenic atom found, not an organometalloidal compound"
    
    # Check for organic framework (at least one carbon is attached to the arsenic)
    organic_framework = any(
        atom.GetAtomicNum() == 6 and 
        any(neighbor.GetAtomicNum() == 33 for neighbor in atom.GetNeighbors())
        for atom in mol.GetAtoms()
    )
    if organic_framework:
        return True, "Arsenic bonded to an organic structure, classified as organometalloidal"
    else:
        return False, "Arsenic not adequately bonded to an organic structure, classify with caution"

# Example Usage:
# print(is_organometalloidal_compound("C[As](O)(O)=O"))  # Should return True with a valid reason