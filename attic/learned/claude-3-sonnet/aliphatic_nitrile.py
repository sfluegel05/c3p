"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:22677 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyano group (C≡N) - more specific pattern that ensures proper connectivity
    cyano_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if cyano_pattern is None:
        return False, "Error in SMARTS pattern"
        
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    
    if not cyano_matches:
        return False, "No cyano (C≡N) group found"
    
    # For each cyano group, check if it's attached to an aliphatic carbon
    for match in cyano_matches:
        nitrile_c = mol.GetAtomWithIdx(match[0])  # Get the carbon of C≡N
        
        # Get the carbon atom attached to the nitrile carbon
        for neighbor in nitrile_c.GetNeighbors():
            if neighbor.GetAtomicNum() != 7:  # Skip the nitrogen atom of C≡N
                # Check if this carbon is NOT part of an aromatic system
                if not neighbor.GetIsAromatic():
                    # Check if the parent carbon is not in a ring or is in an aliphatic ring
                    ring_info = mol.GetRingInfo()
                    if not ring_info.IsAtomInRingOfSize(neighbor.GetIdx(), 6) or not neighbor.GetIsAromatic():
                        # Additional check to ensure the whole connected system is not aromatic
                        aromatic_system = False
                        for next_neighbor in neighbor.GetNeighbors():
                            if next_neighbor.GetIsAromatic():
                                aromatic_system = True
                                break
                        
                        if not aromatic_system:
                            return True, "Contains cyano group (C≡N) attached to aliphatic carbon"
    
    return False, "No cyano groups attached to aliphatic carbons found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:22677',
        'name': 'aliphatic nitrile',
        'definition': 'Any nitrile derived from an aliphatic compound.',
        'parents': ['CHEBI:33577']
    }
}