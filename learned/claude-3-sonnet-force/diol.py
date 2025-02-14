"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:16198 diol
A compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
Aliphatic diols are also called glycols.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    n_hydroxy = len(hydroxy_matches)

    if n_hydroxy == 2:
        return True, "Contains exactly 2 hydroxy groups"
    elif n_hydroxy < 2:
        return False, f"Contains only {n_hydroxy} hydroxy groups, need 2"
    else:
        return False, f"Contains more than 2 hydroxy groups ({n_hydroxy})"

    # Additional checks could include:
    # - Check molecular weight range
    # - Check for specific functional groups or substructures
    # - Check for other properties like logP, rotatable bonds, etc.