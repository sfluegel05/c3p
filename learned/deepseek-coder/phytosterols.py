"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26195 phytosterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol, occurring in plants, with variations in carbon side chains and/or presence/absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general steroid nucleus pattern (tetracyclic ring system with a hydroxyl group at the 3-position)
    steroid_nucleus_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C")
    if steroid_nucleus_pattern is None:
        return False, "Failed to create steroid nucleus pattern"
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "No steroid nucleus found (tetracyclic ring system with a hydroxyl group at the 3-position)"

    # Check for a side chain at the 17-position (typical of sterols)
    # The side chain can vary in length and degree of unsaturation
    side_chain_pattern = Chem.MolFromSmarts("[C@@]12CC[C@@H]3[C@@]4(CC[C@H](C4)O)C[C@@H]1CC[C@@]2(C)CC")
    if side_chain_pattern is None:
        return False, "Failed to create side chain pattern"
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No side chain found at the 17-position"

    # Check molecular weight (phytosterols typically have a higher molecular weight than cholesterol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a typical phytosterol"

    return True, "Contains a steroid nucleus with a hydroxyl group at the 3-position and a side chain at the 17-position"