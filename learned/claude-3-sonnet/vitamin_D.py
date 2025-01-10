"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for characteristic seco-steroid core with triene system
    # Pattern matches the characteristic vitamin D triene system and ring structure
    vitamin_d_pattern = Chem.MolFromSmarts("[CH2,CH]-,=C-,=C-,=C-[C,c]~[C,c]~[C,c]~[C,c]")
    if not mol.HasSubstructMatch(vitamin_d_pattern):
        return False, "Missing characteristic vitamin D triene system"

    # Count hydroxyl groups - vitamin D compounds typically have at least one
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 1:
        return False, "No hydroxyl groups found"

    # Check for basic carbon skeleton size (vitamin D should have significant size)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Carbon skeleton too small for vitamin D"

    # Check for characteristic cyclohexane ring with specific substitution pattern
    cyclohexane_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "Missing characteristic ring structure"

    # Check molecular weight - vitamin D compounds typically between 350-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 750:  # Extended range to account for derivatives
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for vitamin D"

    # Check for characteristic branching pattern
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][CH]")
    if len(mol.GetSubstructMatches(methyl_branch_pattern)) < 2:
        return False, "Missing characteristic methyl branching"

    # Additional check for conjugated double bond system
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated double bond system"

    # If all checks pass, it's likely a vitamin D compound
    reason = ("Contains characteristic vitamin D features: seco-steroid core, "
              f"hydroxyl groups ({hydroxyl_matches}), triene system, and appropriate size")
    return True, reason