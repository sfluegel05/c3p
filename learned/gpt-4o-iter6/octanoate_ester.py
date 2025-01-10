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
        return False, "Invalid SMILES string"
    
    # Updated SMARTS patterns for octanoyl ester
    # Matching exact 8-carbon chain variations with more general branching inclusion
    octanoate_patterns = [
        Chem.MolFromSmarts("C(=O)O[C][C][C][C][C][C][C][C]"),  # Exact generic 8-carbon chain
    ]
    
    for pattern in octanoate_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            for match in matches:
                # Ensure the matched substructure is correctly 8-c chained and checks:
                carbons = [atom for atom in match if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6]
                if len(carbons) == 8:
                    return True, f"Contains octanoyl ester moiety at positions: {matches}"
    
    return False, "No valid octanoate ester linkage with an 8-carbon chain found"