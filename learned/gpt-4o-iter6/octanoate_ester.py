"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester includes an octanoyl group, which is a characteristic 
    C(=O)O group linked to an 8-carbon alkyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # SMARTS pattern for octanoyl ester: [C](=O)OCCCCCCCC
    octanoate_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCC")
    
    if not octanoate_pattern:
        return (None, "Invalid octanoate SMARTS pattern")
    
    matches = mol.GetSubstructMatches(octanoate_pattern)
    
    if matches:
        return True, f"Contains octanoyl ester moiety at positions: {matches}"
    
    return False, "No valid octanoate ester linkage with an 8-carbon chain found"