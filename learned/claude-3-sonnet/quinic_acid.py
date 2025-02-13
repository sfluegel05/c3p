"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:16857 quinic acid
A cyclitol carboxylic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid, with a cyclohexane core,
    at least one hydroxyl group, and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclohexane core
    cyclohexane_pattern = Chem.MolFromSmarts("[C&R1]1[C&R1][C&R1][C&R1][C&R1][C&R1]1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane core found"

    # Look for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Count carbons, oxygens, and hydrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

    if c_count < 6 or c_count > 12:
        return False, "Carbon count outside expected range for quinic acid"
    if o_count < 4 or o_count > 8:
        return False, "Oxygen count outside expected range for quinic acid"
    if h_count < 8 or h_count > 14:
        return False, "Hydrogen count outside expected range for quinic acid"

    # Check molecular weight - quinic acids typically <300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, "Molecular weight too high for quinic acid"

    return True, "Contains cyclohexane core with at least one hydroxyl group and a carboxylic acid"