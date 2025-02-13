"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:50047 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is any saponin derived from a hydroxysteroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclopentanoperhydrophenanthrene scaffold (core steroid structure)
    steroid_core_pattern = Chem.MolFromSmarts("[C&r5,r6]1[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6]1[C&r6]2[C&r5,r6][C&r5,r6][C&r5,r6][C&r5,r6][C&r6]2")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"
    
    # Check for an additional fused ring (5, 6, or 7 carbons)
    ring_sizes = [5, 6, 7]
    has_additional_ring = any(mol.HasSubstructMatch(Chem.MolFromSmarts(f"[C&r{size}]")) for size in ring_sizes)
    if not has_additional_ring:
        return False, "No additional fused ring found"
    
    # Check for OH groups (hydroxysteroid)
    n_oh = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(bond.GetBondTypeAsDouble() < 2 for bond in atom.GetBonds()) > 0)
    if n_oh < 1:
        return False, "No hydroxyl groups found (not a hydroxysteroid)"
    
    # Check for glycosidic bonds (saponin)
    glycosidic_pattern = Chem.MolFromSmarts("[O-;!$(*=,#[!$(*=-,#[!$(*=-,-,-)])>1*]=C")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic bonds found (not a saponin)"
    
    # Typical molecular weight range for steroid saponins
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 2000:
        return False, f"Molecular weight {mol_wt:.0f} Da outside typical range for steroid saponins"
    
    return True, "Contains steroid backbone with hydroxyl groups and glycosidic bonds"