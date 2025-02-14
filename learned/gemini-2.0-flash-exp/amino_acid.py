"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.
    This function explicitly looks for a carbon connecting the two groups, or a simple connection
    through a small chain or ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         tuple(bool, str): True if molecule is an amino acid, False otherwise and the reason.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified SMARTS for amino group (any N with 3 bonds) and carboxyl (carboxylic acid and carboxylate)
    amino_pattern = Chem.MolFromSmarts("[NX3]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX1-,OX2H0,OX2H1]")

    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"


    # Check if the amino and carboxylic groups are connected by a maximum of 5 bonds
    connected_pattern = Chem.MolFromSmarts("[NX3]~[#6]1~[#6]~[#6]~[#6]~[CX3](=[OX1,OX2])[OX1-,OX2H0,OX2H1]")
    connected_pattern2 = Chem.MolFromSmarts("[NX3]~[#6]~[#6]~[#6]~[CX3](=[OX1,OX2])[OX1-,OX2H0,OX2H1]")
    connected_pattern3 = Chem.MolFromSmarts("[NX3]~[#6]~[#6]~[CX3](=[OX1,OX2])[OX1-,OX2H0,OX2H1]")
    connected_pattern4 = Chem.MolFromSmarts("[NX3]~[#6]~[CX3](=[OX1,OX2])[OX1-,OX2H0,OX2H1]")
    connected_pattern5 = Chem.MolFromSmarts("[NX3]~[CX3](=[OX1,OX2])[OX1-,OX2H0,OX2H1]")

    if mol.HasSubstructMatch(connected_pattern) or \
       mol.HasSubstructMatch(connected_pattern2) or \
       mol.HasSubstructMatch(connected_pattern3) or \
       mol.HasSubstructMatch(connected_pattern4) or \
       mol.HasSubstructMatch(connected_pattern5):
         return True, "Contains both amino and carboxylic acid groups with a simple connection"

    return False, "Contains amino and carboxylic acid groups but not in an amino acid-like structure"