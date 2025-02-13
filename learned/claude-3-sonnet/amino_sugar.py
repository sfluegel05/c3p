"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:36973 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino groups (-NH2 or -NHR)
    has_amino = any(atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() >= 1 for atom in mol.GetAtoms())
    if not has_amino:
        return False, "No amino groups found"
    
    # Look for sugar backbone (ring with multiple -OH groups)
    sugar_pattern = Chem.MolFromSmarts("[OX2][CX4]1[CX4]([CX4]([OX2])[CX4]([OX2])[CX4]1)[CX4]([OX2])[CX4]([OX2])")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar backbone found"
    
    # Count hydroxy and amino groups
    n_oh = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    n_nh = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() >= 1)
    
    # Amino sugars should have at least 1 amino and multiple hydroxy groups
    if n_nh < 1 or n_oh < 2:
        return False, "Insufficient amino and/or hydroxy groups"
    
    return True, "Contains a sugar backbone with one or more hydroxy groups replaced by amino groups"