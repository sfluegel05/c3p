"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: hydroxy fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid (long aliphatic chain with a terminal carboxylic acid group)
    carrying one or more hydroxy (-OH) substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) == 0:
        return False, "No carboxylic acid group found"
    
    # Look for hydroxy groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Exclude the hydroxy group in the carboxylic acid
    carboxylic_acid_oh_atoms = [match[2] for match in carboxylic_acid_matches]
    hydroxy_oh_atoms = [match[0] for match in hydroxy_matches]
    # Hydroxy groups excluding those in carboxylic acids
    additional_hydroxy_atoms = [atom_idx for atom_idx in hydroxy_oh_atoms if atom_idx not in carboxylic_acid_oh_atoms]
    
    if len(additional_hydroxy_atoms) == 0:
        return False, "No hydroxy groups found outside carboxylic acid"
    
    # Optionally, check for long aliphatic chain (e.g., at least 4 carbons)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Too few carbon atoms ({num_carbons}) for a fatty acid"
    
    # Ensure that hydroxy groups are attached to carbon chain
    # For each additional hydroxy group, check that it's connected to a carbon atom
    for oh_idx in additional_hydroxy_atoms:
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        neighbors = oh_atom.GetNeighbors()
        if not any(nb.GetAtomicNum() == 6 for nb in neighbors):
            return False, "Hydroxy group not attached to carbon atom"
    
    return True, "Molecule is a hydroxy fatty acid with at least one hydroxy group"

__metadata__ = {
    'chemical_class': {
        'name': 'hydroxy fatty acid',
        'definition': 'Any fatty acid carrying one or more hydroxy substituents.'
    }
}