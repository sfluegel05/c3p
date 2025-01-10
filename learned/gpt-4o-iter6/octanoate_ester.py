"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester includes an octanoyl group, a C(=O)O linked to an 8-carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool, str: True and the reason if molecule is an octanoate ester; otherwise False and reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # SMARTS patterns for octanoyl ester - allowing branching: could be flexible
    octanoate_patterns = [
        Chem.MolFromSmarts("C(=O)OCCCCCCCC"),  # Linear chain
        Chem.MolFromSmarts("C(=O)OCCCCCC(C)C"),  # Branching at end
        Chem.MolFromSmarts("C(=O)O[CH2]CCCCC[CH2]C"),  # Ditto, chain can vary
        Chem.MolFromSmarts("C(=O)OCC(C)CCCC"),  # Midchain branch
    ]
    
    for pattern in octanoate_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            return True, f"Contains octanoyl ester moiety at positions: {matches}"
    
    return False, "No valid octanoate ester linkage with an 8-carbon chain found"