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
    partially or completely hydrogenated, with ring expansions, contractions, and bond scissions
    allowed. Methyl groups at C-10 and C-13 are common, and an alkyl group at C-17 is often present.

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
    
    # Look for cyclopenta[a]phenanthrene backbone (more flexible pattern)
    backbone_pattern = Chem.MolFromSmarts("[C&r5,r6]1[C&r5,r6]2[C&r5,r6]3[C&r5,r6]4[C&r5,r6]5[C&r5,r6]6[C&r5,r6]7[C&r5,r6]8[C&r5,r6]9[C&r5,r6]%91%81%71%61%51%41%31%21%11")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No cyclopenta[a]phenanthrene backbone found"
    
    # Consider additional features
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    cyclopentane_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) == 5)
    cyclohexane_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) == 6)
    
    if mol_weight > 200 and hydroxyl_groups > 0 and (cyclopentane_rings > 0 or cyclohexane_rings > 0):
        return True, "Contains cyclopenta[a]phenanthrene backbone (with possible ring expansions, contractions, or bond scissions), molecular weight > 200 Da, at least one hydroxyl group, and at least one cyclopentane or cyclohexane ring"
    else:
        return False, "Does not meet the criteria for a steroid"