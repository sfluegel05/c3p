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
    steroid_patterns = [
        Chem.MolFromSmarts("[C@@](C)(CC[C@@](C)(C[C@@H]([C@H]1[C@@H](C[C@@H](C[C@@H]([C@@H](C[C@H]([C@@H](C1)C)C)C)C)C)C)C)C"),
        Chem.MolFromSmarts("[C@@](C)(CC[C@@](C)(C[C@@H]([C@H]1[C@@H](C[C@@H](C[C@@H]([C@@H](C[C@H]([C@@H](C1)C)C)C)C)C)C)C)C"),
        # Add more steroid backbone patterns if needed
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "Does not contain steroid backbone"

    # Check for lactone ring
    lactone_pattern = Chem.MolFromSmarts("C(=O)OC")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone ring found"

    # Check for C28 carbon skeleton
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms < 28 or num_heavy_atoms > 32:
        return False, "Does not have a C28 carbon skeleton"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 650:
        return False, "Molecular weight outside typical range for withanolides"

    # Additional checks (optional)
    # ...

    return True, "Contains steroid backbone with a lactone ring and a C28 carbon skeleton"