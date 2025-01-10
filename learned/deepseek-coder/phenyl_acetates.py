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

    # Check if the acetate ester is directly attached to the phenyl ring
    # We look for a pattern where the oxygen of the acetate ester is directly attached to a carbon in the phenyl ring
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1[OX2][CX3](=[OX1])[CX4H3]")
    if not mol.HasSubstructMatch(phenyl_acetate_pattern):
        # Check if the acetate ester is attached to a substituted phenyl ring
        substituted_phenyl_acetate_pattern = Chem.MolFromSmarts("c1c([*])cccc1[OX2][CX3](=[OX1])[CX4H3]")
        if not mol.HasSubstructMatch(substituted_phenyl_acetate_pattern):
            return False, "Acetate ester group not attached to phenyl ring"

    # Ensure that the acetate ester is not part of a larger ester or complex structure
    # We check that the acetate ester is not attached to another ester or complex group
    complex_ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4H3].[OX2][CX3](=[OX1])[CX4H3]")
    if mol.HasSubstructMatch(complex_ester_pattern):
        # Check if the complex ester is part of a larger structure
        # If the molecule has more than one acetate ester, it should still be considered a phenyl acetate if at least one is directly attached to a phenyl ring
        # So we need to ensure that at least one acetate ester is directly attached to a phenyl ring
        # We will count the number of acetate esters directly attached to phenyl rings
        phenyl_acetate_matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
        substituted_phenyl_acetate_matches = mol.GetSubstructMatches(substituted_phenyl_acetate_pattern)
        if len(phenyl_acetate_matches) + len(substituted_phenyl_acetate_matches) == 0:
            return False, "Acetate ester is part of a larger ester or complex structure"

    return True, "Contains a phenyl group with an acetate ester group attached"