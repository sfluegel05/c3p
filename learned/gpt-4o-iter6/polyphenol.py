"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is defined as having 2 or more benzene rings each with at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define pattern for benzene ring (aromatic ring)
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    
    # Identify benzene rings in the molecule
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    # Must have at least 2 benzene rings
    if len(benzene_matches) < 2:
        return False, "Fewer than 2 benzene rings found"

    # Check each benzene ring for at least one hydroxy substitution
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    benzene_with_hydroxy_count = 0
    for match in benzene_matches:
        # Check if any atom in benzene match has a hydroxy neighbor
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.HasSubstructMatch(hydroxy_pattern):
                    benzene_with_hydroxy_count += 1
                    break
    
    # Each benzene ring must have at least one hydroxy substitution
    if benzene_with_hydroxy_count < len(benzene_matches):
        return False, "Not all benzene rings have a hydroxy substitution"

    return True, "Contains 2 or more benzene rings each with at least one hydroxy group"