"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is characterized by having an acetate ester linkage 
    attached to a phenol group or phenyl group.

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

    # Identify acetate ester group (-O-C(=O)CH3)
    acetate_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    acetate_matches = mol.GetSubstructMatches(acetate_ester_pattern)
    if not acetate_matches:
        return False, "No acetate ester linkage found"

    # Identify phenol/phenyl group pattern (aromatic ring with -OH)
    phenyl_esters_pattern = Chem.MolFromSmarts("c1ccc(cc1)OC(=O)C")
    if mol.HasSubstructMatch(phenyl_esters_pattern):
        return True, "Contains acetate ester linkage attached to a phenol/phenyl group"

    return False, "Ester linkage not appropriately attached to a phenol/phenyl group"