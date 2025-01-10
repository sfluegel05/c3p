"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more specific substructure patterns for bisbenzylisoquinoline alkaloids
    ether_bridge_pattern = Chem.MolFromSmarts("COc1ccc(O)cc1")  # A more complex ether linkage
    isoquinoline_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc-2c1")  # Isoquinoline ring with nitrogen
    
    # Checking for two aromatic ether linkages
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_matches) < 2:
        return False, f"Found {len(ether_matches)} ether bridges, need at least 2"
    
    # Checking for isoquinoline units
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(isoquinoline_matches) < 2:
        return False, f"Found {len(isoquinoline_matches)} isoquinoline units, need at least 2"

    # Ensure both components are linked through additional ether bridges
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("c1(COc2ccccc2)c(O)c2[nH]ccc-2c1")):
        return False, "No complex ether linked bisbenzylisoquinoline structure detected"

    return True, "Contains characteristic bisbenzylisoquinoline alkaloid structure"