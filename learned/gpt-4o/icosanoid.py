"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are characterized by a C20 backbone with oxidation at various sites.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a C20 carbon backbone
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Expected at least 20 carbons, found {c_count}"
    
    # Look for oxidation patterns: hydroxyls, ketones, epoxy groups
    has_oxidation = False
    
    # Hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    if mol.HasSubstructMatch(hydroxyl_pattern):
        has_oxidation = True
    
    # Keto groups
    keto_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')
    if mol.HasSubstructMatch(keto_pattern):
        has_oxidation = True

    # Epoxide groups
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    if mol.HasSubstructMatch(epoxide_pattern):
        has_oxidation = True
        
    if not has_oxidation:
        return False, "No typical oxidation patterns (hydroxy, keto, epoxy) found"
    
    # Check for unsaturated bonds typical of icosanoids
    cc_double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bond_matches = len(mol.GetSubstructMatches(cc_double_bond_pattern))
    if double_bond_matches < 3:
        return False, f"Expected at least 3 double bonds, found {double_bond_matches}"
    
    return True, "Contains a C20 backbone with oxidation and multiple double bonds, characteristic of icosanoids"