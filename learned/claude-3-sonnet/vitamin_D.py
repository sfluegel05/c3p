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

    # Characteristic seco-steroid core pattern with broken B-ring and triene system
    # This pattern matches the core structure including the A-ring and C/D ring system
    core_pattern = Chem.MolFromSmarts("""
        [CH2,CH]1[CH2,CH][CH2,CH][C,c]2[C,c]1[C,c][C,c][C,c]2
    """)
    
    # Triene system pattern (more flexible to account for variations)
    triene_pattern = Chem.MolFromSmarts("""
        [C,c]=,:[C,c]-,=[C,c]-,=[C,c]
    """)
    
    # Check for core structure or triene system (some flexibility in matching)
    if not (mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(triene_pattern)):
        return False, "Missing characteristic seco-steroid core or triene system"

    # Look for hydroxyl groups - vitamin D compounds typically have at least one
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 1:
        return False, "No hydroxyl groups found"

    # Check molecular weight - expanded range to account for derivatives
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for vitamin D"

    # Count carbons and check for reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:  # Allowing for nor-derivatives
        return False, "Carbon skeleton too small for vitamin D"

    # Check for characteristic side chain pattern
    side_chain_pattern = Chem.MolFromSmarts("""
        [CH2,CH][CH2,CH][CH2,CH][CH,C]([CH3,OH,*])[CH3,OH,*]
    """)
    
    # More specific pattern for the A-ring with potential hydroxyl
    a_ring_pattern = Chem.MolFromSmarts("""
        [CH2,CH]1[CH2,CH][C,c](=,:[C,c])[C,c][CH2,CH]1
    """)

    # Count matching features
    features = 0
    if mol.HasSubstructMatch(side_chain_pattern):
        features += 1
    if mol.HasSubstructMatch(a_ring_pattern):
        features += 1
    if hydroxyl_matches >= 2:  # Common to have multiple hydroxyls
        features += 1
    if mol.HasSubstructMatch(triene_pattern):
        features += 1

    # Require at least 3 characteristic features
    if features < 3:
        return False, "Insufficient vitamin D structural features"

    # Calculate ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring systems for vitamin D"

    # If all checks pass, it's likely a vitamin D compound
    reason = ("Contains characteristic vitamin D features: "
             f"seco-steroid core, {hydroxyl_matches} hydroxyl groups, "
             f"{ring_count} rings, and appropriate structural elements")
    return True, reason