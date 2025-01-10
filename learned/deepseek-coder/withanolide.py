"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:17855 withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a C28 steroid lactone with a modified side chain forming a lactone ring.

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

    # Check molecular formula (C28)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if not formula.startswith('C28'):
        return False, f"Not C28 steroid (formula: {formula})"

    # Check for steroid backbone (4 fused rings)
    # Simplified pattern for steroid core
    steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CC[C@H](C4)[C@@H]3CC[C@@H]2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for lactone ring in side chain
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found in side chain"

    # Check for typical withanolide features
    # 1. Oxygen at C-22 (common in withanolides)
    # 2. Double bond in side chain
    # 3. Oxygen at C-26 (lactone oxygen)
    withanolide_pattern = Chem.MolFromSmarts("[C@@H]1CC[C@H]2[C@@H]3CC[C@H]4CC[C@H](C4)[C@@H]3CC[C@@H]2C1.[CX3](=[OX1])[OX2][CX4]")
    if not mol.HasSubstructMatch(withanolide_pattern):
        return False, "Missing typical withanolide features"

    # Additional check for molecular weight (typical range for withanolides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} outside typical withanolide range"

    return True, "Contains steroid backbone with lactone ring in side chain and typical withanolide features"