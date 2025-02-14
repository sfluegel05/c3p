"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:37707 secondary ammonium ion
An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion contains an NH2+ group attached to two carbon atoms and no other heteroatoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for secondary ammonium ion pattern
    ammonium_pattern = Chem.MolFromSmarts("[NH2+][CX4,CX3]([#6])[CX4,CX3]([#6])[#6]")
    ammonium_matches = mol.GetSubstructMatches(ammonium_pattern)
    
    if not ammonium_matches:
        return False, "No secondary ammonium ion pattern found"
    
    # Check if NH2+ is not attached to heteroatoms other than carbon
    for match in ammonium_matches:
        nitrogen = mol.GetAtomWithIdx(match)
        neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in nitrogen.GetNeighbors()]
        hetero_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() != 6]
        
        if hetero_neighbors:
            continue  # Skip this match if nitrogen is attached to heteroatoms
        
        # Check if NH2+ is not part of a ring
        if nitrogen.IsInRing():
            continue  # Skip this match if nitrogen is part of a ring
        
        # Check if nitrogen has a formal charge of +1
        if nitrogen.GetFormalCharge() != 1:
            continue  # Skip this match if nitrogen is not protonated
        
        return True, "Contains a secondary ammonium ion (NH2+ attached to two carbon atoms and no other heteroatoms)"
    
    return False, "No valid secondary ammonium ion found in the molecule"