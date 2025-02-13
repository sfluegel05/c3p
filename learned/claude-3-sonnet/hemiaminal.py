"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:35841 hemiaminal
A hemiaminal is any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find amino groups
    amino_pattern = Chem.MolFromSmarts("[NX3]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Find hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check if any amino and hydroxy groups are attached to the same carbon
    for amino_idx in amino_matches:
        for hydroxy_idx in hydroxy_matches:
            amino_atom = mol.GetAtomWithIdx(amino_idx)
            hydroxy_atom = mol.GetAtomWithIdx(hydroxy_idx)
            
            if list(amino_atom.GetNeighbors()) & list(hydroxy_atom.GetNeighbors()):
                return True, "Contains amino and hydroxy groups attached to the same carbon"
    
    return False, "No amino and hydroxy groups found attached to the same carbon"