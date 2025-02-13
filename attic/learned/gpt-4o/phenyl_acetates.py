"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is characterized by having an acetate ester linkage 
    attached to a phenol or phenyl group.

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

    # Identify aromatic phenyl ring pattern, tolerance for substitution
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")  # Basic phenyl ring
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)
    
    if not phenyl_matches:
        return False, "No phenyl group found"

    # Identify acetate ester group (-O-C(=O)-C)
    acetate_ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    acetate_matches = mol.GetSubstructMatches(acetate_ester_pattern)
    
    if not acetate_matches:
        return False, "No acetate ester linkage found"

    # Verify that an acetate ester is attached directly to the phenyl group
    for phenyl_match in phenyl_matches:
        for idx in phenyl_match:  # each atom in the phenyl group
            atom = mol.GetAtomWithIdx(idx)
            # Check neighboring atoms to see if any match the start of the acetate
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if any(neighbor_idx in acetate_match for acetate_match in acetate_matches):
                    return True, "Contains acetate ester linkage attached to a phenyl group"
    
    return False, "Ester linkage not attached to phenyl group"