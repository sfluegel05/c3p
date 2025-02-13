"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol contains 2 or more benzene rings, each substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for benzene ring and hydroxy group
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')  # pattern for benzene ring
    hydroxy_pattern = Chem.MolFromSmarts('c[OH]')    # Phenolic -OH group attached to an aromatic carbon
    
    # Find all benzene rings in the molecule
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if len(benzene_matches) < 2:
        return False, "Fewer than two benzene rings"
    
    # Check if each benzene ring has at least one hydroxy group
    rings_with_hydroxy = 0
    for benzene_ring in benzene_matches:
        # Create a substructure for the ring
        ring_atoms = set(benzene_ring)
        ring_molecule = Chem.PathToSubmol(mol, benzene_ring)
        
        # Find all hydroxy groups in the sub molecule
        hydroxy_matches = ring_molecule.GetSubstructMatches(hydroxy_pattern)
        if hydroxy_matches:
            rings_with_hydroxy += 1
    
    if rings_with_hydroxy < 2:
        return False, f"Only {rings_with_hydroxy} benzene rings have a hydroxy group"
        
    return True, "Contains 2 or more benzene rings each with an OH group"