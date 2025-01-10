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

    # Look for cyano group (C≡N)
    cyano_pattern = Chem.MolFromSmarts("[C#N]")
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    
    if not cyano_matches:
        return False, "No cyano (C≡N) group found"
    
    # For each cyano group, check if it's attached to an aliphatic carbon
    for match in cyano_matches:
        nitrile_n = mol.GetAtomWithIdx(match[0])
        
        # Get the carbon atom attached to the nitrile nitrogen
        for neighbor in nitrile_n.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if this carbon is part of an aromatic system
                if neighbor.GetIsAromatic():
                    continue  # Skip this cyano group if it's attached to aromatic carbon
                
                # If we found at least one cyano group attached to non-aromatic carbon,
                # this is an aliphatic nitrile
                return True, "Contains cyano group (C≡N) attached to aliphatic carbon"
    
    return False, "All cyano groups are attached to aromatic carbons or no valid cyano groups found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:22677',
        'name': 'aliphatic nitrile',
        'definition': 'Any nitrile derived from an aliphatic compound.',
        'parents': ['CHEBI:33577']
    }
}