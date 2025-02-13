"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with a characteristic benzopyran-4-one core structure.

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

    # Basic benzopyran-4-one core structure of flavonoids
    benzopyrone_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc21")
    flavonoid_core = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc21")
    
    # More specific patterns for different types of flavonoids
    flavone_pattern = Chem.MolFromSmarts("O=C1CC(=C(c2ccccc2))[O]c2ccccc21")
    
    # Check for presence of core structure
    if not (mol.HasSubstructMatch(benzopyrone_pattern) or 
           mol.HasSubstructMatch(flavonoid_core) or
           mol.HasSubstructMatch(flavone_pattern)):
        return False, "Missing benzopyran-4-one core structure"

    # Count oxygen atoms (should have multiple due to OH/OMe groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for anthoxanthin"

    # Look for hydroxyl groups (common in anthoxanthins)
    oh_pattern = Chem.MolFromSmarts("[OH]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Look for methoxy groups (common in anthoxanthins)
    ome_pattern = Chem.MolFromSmarts("[OH0]C")
    ome_matches = mol.GetSubstructMatches(ome_pattern)
    
    total_substituents = len(oh_matches) + len(ome_matches)
    if total_substituents < 1:
        return False, "Missing characteristic OH/OMe substituents"

    # Look for aromatic rings (should have at least 2)
    aromatic_rings = 0
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic_rings += 1
    if aromatic_rings < 8:  # Each benzene ring contributes 6 aromatic atoms
        return False, "Insufficient aromatic character"

    # Check for ketone group
    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    # Additional check for carbon count (flavonoids typically have at least 15 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Insufficient carbon atoms for flavonoid structure"

    return True, "Contains benzopyran-4-one core with appropriate substituents and aromatic character"