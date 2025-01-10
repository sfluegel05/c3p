"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    Bisbenzylisoquinoline alkaloids are characterized by complex polycyclic structures with multiple
    benzylisoquinoline-like units linked often via ether bonds.

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
    
    # Define a general benzylisoquinoline-like pattern (aromatic rings with N heterocycle)
    benzylisoquinoline_pattern = Chem.MolFromSmarts("Nc1c(cc(OC)c(OC)c1)CC2")
    
    # Check for benzylisoquinoline-like patterns at least twice
    benzylisoquinoline_matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(benzylisoquinoline_matches) < 2:
        return False, f"Found {len(benzylisoquinoline_matches)} benzylisoquinoline-like units, need at least 2"

    # Define smarter ether bond linkers connectivity, not requiring fixed bridge counts
    # Just ensure presence between rings
    ether_pattern = Chem.MolFromSmarts("cOc") 
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 1:
        return False, f"Found {len(ether_matches)} ether bridges, need at least 1"

    return True, "Contains characteristic bisbenzylisoquinoline alkaloid structure"