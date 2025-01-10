"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    Bisbenzylisoquinoline alkaloids are characterized by two benzylisoquinoline units 
    linked by ether bonds with complex aromatic ring systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for benzylisoquinoline unit
    benzylisoquinoline_pattern = Chem.MolFromSmarts("c1ccc2Nc3c(ccc(O)c3)Cc2c1")
    
    # Define SMARTS for ether bridge (between aromatic rings)
    ether_bridge_pattern = Chem.MolFromSmarts("c1cc(OC)cc1")
    
    # Check for two benzylisoquinoline patterns
    benzylisoquinoline_matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(benzylisoquinoline_matches) < 2:
        return False, f"Found {len(benzylisoquinoline_matches)} benzylisoquinoline units, need at least 2"
    
    # Check for ether bridges (requires two such bridges)
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_matches) < 2:
        return False, f"Found {len(ether_matches)} ether bridges, need at least 2"
    
    return True, "Contains characteristic bisbenzylisoquinoline alkaloid structure"