"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:51800 withanolide

A withanolide is any steroid lactone that is a C28 steroid with a modified side chain
forming a lactone ring and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a withanolide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@@](C)(CC[C@@](C)(C[C@@H]([C@H]1[C@@H](C[C@@H](C[C@@H]([C@@H](C[C@H]([C@@H](C1)C)C)C)C)C)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid backbone"

    # Check for lactone ring
    lactone_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for C28 carbon skeleton
    if mol.GetNumHeavyAtoms() != 28:
        return False, "Does not have a C28 carbon skeleton"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula != "C28H38O4" and formula != "C28H40O4":
        return False, "Molecular formula not consistent with withanolides"

    # Additional checks (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 600:
        return False, "Molecular weight outside typical range for withanolides"

    return True, "Contains steroid backbone with a lactone ring and a C28 carbon skeleton"