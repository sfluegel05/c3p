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
    Anthoxanthins are flavonoid pigments that include flavones and flavonols.

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

    # More flexible core patterns for flavonoids
    # Basic chromone (benzopyran-4-one) core
    chromone_pattern = Chem.MolFromSmarts("O=C1CCOc2ccccc12")
    
    # Flavone core (more specific)
    flavone_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    
    # Flavonol core (with 3-OH group)
    flavonol_pattern = Chem.MolFromSmarts("O=C1C(O)C(c2ccccc2)Oc2ccccc12")
    
    # Alternative flavonoid core pattern
    alt_flavonoid = Chem.MolFromSmarts("O=C1CC(=C)Oc2ccccc12")

    # Check for presence of any core structure
    has_core = any([
        mol.HasSubstructMatch(pattern) 
        for pattern in [chromone_pattern, flavone_pattern, flavonol_pattern, alt_flavonoid]
        if pattern is not None
    ])
    
    if not has_core:
        return False, "Missing flavonoid core structure"

    # Count oxygen atoms (should have multiple due to OH/OMe groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient oxygen atoms for anthoxanthin"

    # Look for hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OH]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Look for methoxy groups
    ome_pattern = Chem.MolFromSmarts("OC")
    ome_matches = mol.GetSubstructMatches(ome_pattern)
    
    # Look for glycoside patterns (common in anthoxanthins)
    glycoside_pattern = Chem.MolFromSmarts("OC1OCC(O)C(O)C1O")
    has_glycoside = mol.HasSubstructMatch(glycoside_pattern) if glycoside_pattern else False
    
    total_substituents = len(oh_matches) + len(ome_matches)
    if total_substituents < 1 and not has_glycoside:
        return False, "Missing characteristic OH/OMe/glycoside substituents"

    # Look for aromatic systems
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 6:  # At least one aromatic ring
        return False, "Insufficient aromatic character"

    # Check for ketone group (part of the chromone system)
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing ketone group"

    # Count carbons (flavonoids typically have at least 15 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:  # Lowered threshold to catch more variants
        return False, "Insufficient carbon atoms for flavonoid structure"

    # Additional check for characteristic double bond
    double_bond = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond):
        return False, "Missing characteristic double bond"

    return True, "Contains flavonoid core structure with appropriate substituents"