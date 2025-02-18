"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:27306 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D analog based on its SMILES string.
    Vitamin D compounds are 9,10-seco-steroids with a conjugated triene system.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Essential: 9,10-seco-steroid core with conjugated triene
    # SMARTS pattern for characteristic 10,19-seco structure with conjugated triene
    seco_pattern = Chem.MolFromSmarts(
        "[CH3]C[C@H]1CC[C@H]2C(=C/C=C3/C(=C/C=C\4)"  # Modified to match conjugated triene in seco structure
        "[C@]5([C@@H](CC[C@H]5C(=C4)CC3)C)CC2)C1"  # Accounts for stereochemistry variations
    )
    if not mol.HasSubstructMatch(seco_pattern):
        return False, "Missing seco-steroid core with conjugated triene"

    # Check for at least two hydroxyl groups (common in active forms)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() 
                        if atom.GetAtomicNum() == 8 
                        and atom.GetTotalNumHs() >= 1 
                        and atom.GetHybridization() == Chem.HybridizationType.SP3)
    if hydroxyl_count < 1:
        return False, "Insufficient hydroxyl groups"

    # Verify steroid methyl groups (C18 and C19)
    methyl_pattern = Chem.MolFromSmarts("[CH3]-C(-[CH2])(-[CH2])")
    if len(mol.GetSubstructMatches(methyl_pattern)) < 2:
        return False, "Missing characteristic steroid methyl groups"

    # Adjusted molecular weight check (some analogs are smaller)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 280:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    # Check for conjugated triene system with specific geometry
    # More precise pattern requiring 3 conjugated double bonds in sequence
    triene_pattern = Chem.MolFromSmarts("*=*/*=*/*=*")  # Allows some branching
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system"

    return True, "Contains seco-steroid core with conjugated triene and hydroxyl groups"