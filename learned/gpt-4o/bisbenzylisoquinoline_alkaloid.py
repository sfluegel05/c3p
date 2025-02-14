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

    # SMARTS patterns for benzylisoquinoline units. We refine it to match the common structure:
    # Benzylisoquinoline unit often looks like this with variability at positions:
    benzylisoquinoline_pattern = Chem.MolFromSmarts('c1ccc2cc3c(ccc4c3c(cn2)c(c1)c4)') 
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, "Does not contain two benzylisoquinoline units"

    # Refined search for ether bridges between two benzylisoquinoline structures
    # Allow for spacing in the ether bridge in complex structures
    ether_bridge_pattern = Chem.MolFromSmarts('c-O-c')
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_matches) == 0:
        return False, "No ether bridge detected between benzylisoquinoline units"

    # Allow for International Union of Pure and Applied Chemistry (IUPAC) variability; e.g., methylenedioxy bridges, etc.
    methylenedioxy_pattern = Chem.MolFromSmarts('COC1COC1')
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)
    
    # Check for bridges (include flexible pattern with aromatic possibility)
    # This can catch additional setups like those with methylenedioxy and allow for common covalent linkages
    flexible_bridge_pattern = Chem.MolFromSmarts('c-c-c')
    flexible_bridge_matches = mol.GetSubstructMatches(flexible_bridge_pattern)
    
    if len(methylenedioxy_matches) > 0 or len(flexible_bridge_matches) > 0:
        return True, "Contains common bridging motifs for bisbenzylisoquinoline alkaloid"

    return False, "No adequate bridging or secondary bridging units found for bisbenzylisoquinoline alkaloid classification"