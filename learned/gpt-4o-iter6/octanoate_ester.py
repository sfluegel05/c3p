"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester features the presence of the octanoyl group `(C(=O)OCCCCCCCC)`.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for octanoyl ester linkage (-C(=O)OCCCCCCCC)
    octanoyl_ester_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    if not mol.HasSubstructMatch(octanoyl_ester_pattern):
        return False, "No octanoyl ester linkage found"
    
    return True, "Contains octanoyl ester linkage"