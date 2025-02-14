"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:35349 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as any compound based on the cyclopenta[a]phenanthrene carbon skeleton,
    with methyl groups at C-10 and C-13 being preferred, and an alkyl group at C-17 being
    common but not required. Ring expansions, contractions, and bond scissions are allowed.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopenta[a]phenanthrene backbone (more general pattern)
    backbone_pattern = Chem.MolFromSmarts("[C@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)CC2)CC1")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No cyclopenta[a]phenanthrene backbone found"
    
    # Check for methyl groups at C-10 and C-13 (preferred but not required)
    c10_pattern = Chem.MolFromSmarts("[C@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)CC2)C(C)CC1")
    c13_pattern = Chem.MolFromSmarts("[C@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C)CC2)CC(C)C1")
    c10_match = mol.HasSubstructMatch(c10_pattern)
    c13_match = mol.HasSubstructMatch(c13_pattern)
    methyl_groups_present = c10_match and c13_match
    
    # Check for alkyl group at C-17 (common but not required)
    c17_pattern = Chem.MolFromSmarts("[C@]12[C@@H]([C@@H]3[C@@H]([C@@H]([C@@H]4[C@@H]([C@H](CC4)CC3)C[C@@H]([C@H](C)CC)C)CC2)CC1")
    c17_match = mol.HasSubstructMatch(c17_pattern)
    
    # Allow for ring expansions, contractions, and bond scissions
    ring_size_check = all(len(ring) >= 5 and len(ring) <= 7 for ring in mol.GetRingInfo().AtomRings())
    disconnected_substructures = len(Chem.GetMolFrags(mol)) > 1
    
    # Consider additional features
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    cyclopentane_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) == 5)
    cyclohexane_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) == 6)
    
    if methyl_groups_present and c17_match and ring_size_check and not disconnected_substructures and mol_weight > 200 and hydroxyl_groups > 0 and (cyclopentane_rings > 0 or cyclohexane_rings > 0):
        return True, "Contains cyclopenta[a]phenanthrene backbone with methyl groups at C-10 and C-13, alkyl group at C-17, no unusual ring sizes or disconnected substructures, molecular weight > 200 Da, at least one hydroxyl group, and at least one cyclopentane or cyclohexane ring"
    elif methyl_groups_present and ring_size_check and not disconnected_substructures and mol_weight > 200 and hydroxyl_groups > 0 and (cyclopentane_rings > 0 or cyclohexane_rings > 0):
        return True, "Contains cyclopenta[a]phenanthrene backbone with methyl groups at C-10 and C-13, no unusual ring sizes or disconnected substructures, molecular weight > 200 Da, at least one hydroxyl group, and at least one cyclopentane or cyclohexane ring"
    elif ring_size_check and not disconnected_substructures and mol_weight > 200 and hydroxyl_groups > 0 and (cyclopentane_rings > 0 or cyclohexane_rings > 0):
        return True, "Contains cyclopenta[a]phenanthrene backbone with no unusual ring sizes or disconnected substructures, molecular weight > 200 Da, at least one hydroxyl group, and at least one cyclopentane or cyclohexane ring"
    else:
        return False, "Does not meet the criteria for a steroid"