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

    # Characteristic triene system pattern
    triene_pattern = Chem.MolFromSmarts("C=CC=CC=C")
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "Missing characteristic triene system"

    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 1:
        return False, "No hydroxyl groups found"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 650:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for vitamin D"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Most vitamin D compounds have at least 20 carbons
        return False, "Carbon skeleton too small for vitamin D"

    # Ring count - vitamin D compounds typically have 3-4 rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2 or ring_count > 5:
        return False, f"Ring count {ring_count} unusual for vitamin D structure"

    # Check for characteristic cyclohexene A-ring pattern
    a_ring_pattern = Chem.MolFromSmarts("C1CCC(=C)CC1")
    
    # Side chain pattern (more flexible to account for variations)
    side_chain_pattern = Chem.MolFromSmarts("CCCC(C)C")

    # Score features
    features = 0
    if mol.HasSubstructMatch(a_ring_pattern):
        features += 1
    if mol.HasSubstructMatch(side_chain_pattern):
        features += 1
    if hydroxyl_matches >= 2:  # Common to have multiple hydroxyls
        features += 1
    if c_count >= 27:  # Full-length vitamin D typically has 27 carbons
        features += 1

    # Check degree of unsaturation
    double_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if double_bonds >= 3:
        features += 1

    # Require at least 3 characteristic features
    if features < 3:
        return False, "Insufficient vitamin D structural features"

    # If all checks pass
    reason = (f"Contains vitamin D characteristics: triene system, {hydroxyl_matches} hydroxyl groups, "
             f"{ring_count} rings, {c_count} carbons, and appropriate structural elements")
    return True, reason