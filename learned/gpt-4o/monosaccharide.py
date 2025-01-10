"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones, with at least three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Identify the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return (False, "Too few carbon atoms; monosaccharides must have at least 3")
    
    # Identify the number of hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group SMARTS pattern
    hydroxyl_matches = mol.GetSubstructMatches(oh_pattern)
    if len(hydroxyl_matches) < 2:
        return (False, f"Insufficient hydroxyl groups found; needed at least 2, found {len(hydroxyl_matches)}")

    # Check for potential aldehyde or ketone functionalities or presence of hemiketals/hemiacetals
    carbonyl_pattern = Chem.MolFromSmarts("[$([CX3]=[OX1]),$([CX3H1][OX2H])]")  # Carbonyl under acyclic or cyclic conditions
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return (False, "No potential carbonyl or cyclic oxy-functional group detected.")
    
    # Verify molecular weight does not exceed typical monosaccharide weight (~300 Da)
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 300:
        return (False, "Molecular weight too high; suggests multi-unit composition.")

    return (True, "Structure is a monosaccharide with potential cyclic or acyclic carbonyl group and sufficient hydroxyls.")