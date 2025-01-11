"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 skeleton, often with hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the C6-C3-C6 skeleton pattern
    # This pattern matches two aromatic rings connected by a three-atom chain, allowing for variations
    skeleton_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4,CX3,O]-[CX4,CX3,O]-[CX4,CX3,O]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt1 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4,CX3,O]-[CX4,CX3,O]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt2 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4,CX3,O]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt3 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4,CX3,O]-[CX4,CX3,O]-[CX4,CX3,O]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt4 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4,CX3,O]-[CX4,CX3,O]-[c]2[c][c][c][c][c]2")
    skeleton_pattern_alt5 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4,CX3,O]-[c]2[c][c][c][c][c]2")

    if not (mol.HasSubstructMatch(skeleton_pattern) or 
            mol.HasSubstructMatch(skeleton_pattern_alt1) or 
            mol.HasSubstructMatch(skeleton_pattern_alt2) or 
            mol.HasSubstructMatch(skeleton_pattern_alt3) or 
            mol.HasSubstructMatch(skeleton_pattern_alt4) or 
            mol.HasSubstructMatch(skeleton_pattern_alt5)):
        return False, "No C6-C3-C6 skeleton found"

    # Count aromatic rings to ensure there are at least two
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, f"Found {aromatic_rings} aromatic rings, need at least 2"

    # Check for hydroxyl groups attached to the aromatic rings
    hydroxyl_pattern = Chem.MolFromSmarts("[c][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 1"

    # Check molecular weight - flavonoids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for flavonoid"

    return True, "Contains C6-C3-C6 skeleton with at least one hydroxyl group"