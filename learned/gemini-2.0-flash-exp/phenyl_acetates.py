"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate has an acetate group attached to a phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phenyl ring
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl ring found"

    # Check for acetate group connected to oxygen
    acetate_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]-[OX2]")
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)

    if not acetate_matches:
        return False, "No acetate group found"

    for match in acetate_matches:
        # The last oxygen in the match is the one that should be attached to phenyl
        oxygen_index = match[2]
        oxygen_atom = mol.GetAtomWithIdx(oxygen_index)
        
        for neighbor in oxygen_atom.GetNeighbors():
           for phenyl_atom_idx in mol.GetSubstructMatch(phenyl_pattern):
               if neighbor.GetIdx() == phenyl_atom_idx:
                   return True, "Contains a phenyl ring with at least one acetate group attached to a phenol oxygen"

    return False, "Acetate is not attached to a phenyl ring"