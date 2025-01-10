"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate derivative based on its SMILES string.
    A phenyl acetate generally has an acetyl group attached to a phenyl ring.

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
    
    # Phenyl ring pattern (benzene)
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl ring found"
    
    # Acetyl group pattern
    acetyl_pattern = Chem.MolFromSmarts("C(=O)C")
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    
    # Check if any acetyl group has a bond to a phenyl carbon
    for match in acetyl_matches:
        acetyl_atom = [atom for atom in match if mol.GetAtomWithIdx(atom).GetSymbol() == "C" and len(mol.GetAtomWithIdx(atom).GetNeighbors()) > 2]
        for atom_idx in acetyl_atom:
            neighbors = mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            if any(neigh.GetSymbol() == "C" and neigh.GetIsAromatic() for neigh in neighbors):
                return True, "Contains phenyl ring with acetyl group attached"
                
    return False, "No acetyl group attached to phenyl ring"