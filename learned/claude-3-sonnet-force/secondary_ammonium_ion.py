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
    A secondary ammonium ion contains an NH2+ group attached to two carbon atoms.

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

    # Look for NH2+ pattern
    ammonium_pattern = Chem.MolFromSmarts("[NH2+]")
    ammonium_matches = mol.GetSubstructMatches(ammonium_pattern)
    
    if not ammonium_matches:
        return False, "No NH2+ group found"
    
    # Check if NH2+ is attached to two carbon atoms
    for match in ammonium_matches:
        nitrogen = mol.GetAtomWithIdx(match)
        neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in nitrogen.GetNeighbors()]
        carbon_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 6]
        
        if len(carbon_neighbors) == 2:
            return True, "Contains an NH2+ group attached to two carbon atoms (secondary ammonium ion)"
    
    return False, "NH2+ group not attached to two carbon atoms"