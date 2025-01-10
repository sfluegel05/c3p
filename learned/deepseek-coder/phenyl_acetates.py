"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester obtained by formal condensation of the carboxy group of acetic acid with the hydroxy group of any phenol.

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

    # Define the acetate ester pattern: -O-C(=O)-CH3
    acetate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4H3]")
    if not mol.HasSubstructMatch(acetate_pattern):
        return False, "No acetate ester group found"

    # Define the phenyl group pattern: benzene ring
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl group found"

    # Check if the acetate ester is attached to the phenyl group
    # We look for a pattern where the oxygen of the acetate ester is directly attached to a carbon in the phenyl ring
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1[OX2][CX3](=[OX1])[CX4H3]")
    if not mol.HasSubstructMatch(phenyl_acetate_pattern):
        return False, "Acetate ester group not attached to phenyl ring"

    return True, "Contains a phenyl group with an acetate ester group attached"