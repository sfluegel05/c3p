"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: CHEBI:38623 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with a flavone or flavonol backbone,
    typically containing multiple hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavone/flavonol backbone (C6-C3-C6 structure with a carbonyl group)
    # More general pattern to capture variations in the backbone
    flavone_pattern = Chem.MolFromSmarts("[O]=C1C=C2C(=C1)C(=O)C=C2")
    flavone_pattern_general = Chem.MolFromSmarts("[O]=C1C=C2C(=C1)C(=O)C=C2")
    if not mol.HasSubstructMatch(flavone_pattern) and not mol.HasSubstructMatch(flavone_pattern_general):
        return False, "No flavone/flavonol backbone found"

    # Count hydroxyl groups (OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxyl_count < 2:
        return False, f"Found {hydroxyl_count} hydroxyl groups, need at least 2"

    # Check for aromatic rings (should have at least 2 aromatic rings)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, f"Found {aromatic_rings} aromatic rings, need at least 2"

    # Check molecular weight (anthoxanthins typically have MW > 200 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for anthoxanthin"

    return True, "Contains flavone/flavonol backbone with multiple hydroxyl groups"