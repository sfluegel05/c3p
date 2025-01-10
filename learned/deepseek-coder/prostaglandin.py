"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:26347 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is a naturally occurring compound derived from the parent C20 acid, prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a cyclopentane ring with at least one substituent
    cyclopentane_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"

    # Check for the presence of at least one carboxyl group (COOH or ester)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for the presence of at least one hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for the presence of at least one double bond (C=C)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond found"

    # Check molecular weight - prostaglandins typically have a molecular weight around 300-600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 650:
        return False, f"Molecular weight {mol_wt:.2f} is outside the typical range for prostaglandins"

    # Count carbons - prostaglandins typically have 20 carbons, but allow some flexibility
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 25:
        return False, f"Number of carbons {c_count} is outside the typical range for prostaglandins"

    # Check for the presence of a specific prostaglandin-like structure
    prostaglandin_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1.[CX3](=O)[OX2H0].[OX2H].[CX3]=[CX3]")
    if not mol.HasSubstructMatch(prostaglandin_pattern):
        return False, "No prostaglandin-like structure found"

    return True, "Contains cyclopentane ring, carboxyl group, hydroxyl group, and double bond, with appropriate molecular weight and carbon count"