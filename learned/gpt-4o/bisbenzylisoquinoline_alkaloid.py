"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determine if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
   
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the bisbenzylisoquinoline alkaloid pattern, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS pattern for the benzylisoquinoline unit
    benzylisoquinoline_pattern = Chem.MolFromSmarts('c1ccc2cc3c(ccc4c3c(cn2)c(c1)C4)') 
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, "Does not contain two benzylisoquinoline units"

    # Check for an ether bridge linking two benzylisoquinoline structures
    ether_bridge_pattern = Chem.MolFromSmarts('c-O-c')  # Simplified to detect C-O-C linkage
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_matches) == 0:
        return False, "No ether bridge detected between benzylisoquinoline units"

    # Methylenedioxy bridge pattern
    methylenedioxy_pattern = Chem.MolFromSmarts('COC1COC1')
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)
    if len(methylenedioxy_matches) > 0:
        return True, "Contains methylenedioxy bridge indicating bisbenzylisoquinoline alkaloid"

    # Further relaxed search for additional common C-C bridging connections
    cc_bridge_pattern = Chem.MolFromSmarts('c1cc2ccccc2c1-c3ccccc3')
    cc_bridge_matches = mol.GetSubstructMatches(cc_bridge_pattern)
    if len(cc_bridge_matches) > 0:
        return True, "Contains C-C bridge indicating bisbenzylisoquinoline alkaloid"

    return True, "Classified as bisbenzylisoquinoline alkaloid based on structural patterns"