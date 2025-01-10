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
    cyano_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if not mol.HasSubstructMatch(cyano_pattern):
        return False, "No cyano (C≡N) group found"
    
    # Get all cyano group matches
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    
    # Pattern for aromatic carbon attached to cyano group
    aromatic_cyano = Chem.MolFromSmarts("a[CX2]#[NX1]")
    
    # Pattern for conjugated system directly attached to cyano
    conjugated_cyano = Chem.MolFromSmarts("[$(C=C),$(C=N),$(C=O)][CX2]#[NX1]")
    
    for match in cyano_matches:
        nitrile_c = mol.GetAtomWithIdx(match[0])  # Get the carbon of C≡N
        
        # Get the atom attached to the nitrile carbon (excluding the nitrogen)
        neighbors = [atom for atom in nitrile_c.GetNeighbors() 
                    if atom.GetAtomicNum() != 7]
        
        if not neighbors:
            continue
            
        attached_atom = neighbors[0]
        
        # Check if the cyano group is attached to an aromatic system
        if attached_atom.GetIsAromatic():
            continue
            
        # If the molecule has the aromatic-cyano pattern at this position, skip
        if aromatic_cyano and mol.HasSubstructMatch(aromatic_cyano):
            continue
            
        # If the cyano group is directly attached to a conjugated system, skip
        if conjugated_cyano and mol.HasSubstructMatch(conjugated_cyano):
            continue
            
        # If we get here, we have an aliphatic nitrile
        return True, "Contains cyano group (C≡N) attached to an aliphatic carbon"
    
    return False, "No aliphatic nitrile groups found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:22677',
        'name': 'aliphatic nitrile',
        'definition': 'Any nitrile derived from an aliphatic compound.',
        'parents': ['CHEBI:33577']
    }
}