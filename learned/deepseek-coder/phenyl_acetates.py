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

    # Define the phenyl acetate pattern: acetate ester attached to a phenyl ring
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1ccccc1[OX2][CX3](=[OX1])[CX4H3]")
    # Define the substituted phenyl acetate pattern: acetate ester attached to a substituted phenyl ring
    substituted_phenyl_acetate_pattern = Chem.MolFromSmarts("c1c([*])cccc1[OX2][CX3](=[OX1])[CX4H3]")

    # Check if the acetate ester is directly attached to the phenyl ring or a substituted phenyl ring
    if not mol.HasSubstructMatch(phenyl_acetate_pattern) and not mol.HasSubstructMatch(substituted_phenyl_acetate_pattern):
        return False, "Acetate ester group not attached to phenyl ring"

    # If the molecule has more than one acetate ester, ensure at least one is attached to a phenyl ring
    acetate_matches = mol.GetSubstructMatches(acetate_pattern)
    if len(acetate_matches) > 1:
        # Count the number of acetate esters attached to phenyl rings
        phenyl_acetate_matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
        substituted_phenyl_acetate_matches = mol.GetSubstructMatches(substituted_phenyl_acetate_pattern)
        if len(phenyl_acetate_matches) + len(substituted_phenyl_acetate_matches) == 0:
            return False, "Multiple acetate esters found, but none attached to phenyl ring"

    return True, "Contains a phenyl group with an acetate ester group attached"