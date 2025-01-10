"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is characterized by having an acetate ester linkage 
    directly attached to a phenol group or phenyl group.

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

    # Identify acetate ester group (-O-C(=O)CH3) or more flexible
    # Include optional esther bonds connecting benzene derivatives
    acetate_ester_pattern = Chem.MolFromSmarts("C(=O)[OX2H1,R]")
    if not mol.HasSubstructMatch(acetate_ester_pattern):
        return False, "No acetate ester linkage found"

    # Identify phenyl/phenol group with ester (-O-C(=O))
    phenyl_esters_pattern = Chem.MolFromSmarts("c1ccc(cc1)OC(=O)")
    phenol_pattern = Chem.MolFromSmarts("Oc1ccccc1OC(=O)")
    
    # Check for either pattern indicating direct attachment
    if mol.HasSubstructMatch(phenyl_esters_pattern) or mol.HasSubstructMatch(phenol_pattern):
        return True, "Contains acetate ester linkage directly attached to a phenyl or phenol group"

    return False, "Ester linkage not appropriately attached to a phenol/phenyl group"