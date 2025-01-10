"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is defined by having an acetate ester linkage attached to a phenol or phenyl group.

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

    # Identify aromatic ring (phenyl group)
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")  # Basic phenyl ring
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl group found"
    
    # Identify acetate ester group (-O-C(=O)-CH3)
    acetate_ester_pattern = Chem.MolFromSmarts("O=C(O)c")  # Acetate ester linkage
    if not mol.HasSubstructMatch(acetate_ester_pattern):
        return False, "No acetate ester linkage found"
    
    # Verify connection of ester to phenyl group
    subms = mol.GetSubstructMatches(phenyl_pattern)
    for sub in subms:
        phenyl = mol.GetAtomWithIdx(sub[0]).GetNeighbors()
        for atom in phenyl:
            if atom.HasSubstructMatch(acetate_ester_pattern):
                return True, "Contains acetate ester linkage attached to a phenyl group"
            
    return False, "Ester linkage not attached to phenyl group"