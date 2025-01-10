"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are characterized by essential fatty acid oxidation products often having C20 or variant backbones with unsaturation and oxidation sites.

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
    
    # Look for larger carbon backbone to allow variation
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16 or c_count > 26:
        return False, f"Expected 16-26 carbons, found {c_count}"
    
    # Look for oxidation patterns: hydroxyls, ketones, epoxide, lactone groups
    has_oxidation = False
    
    # Hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    if mol.HasSubstructMatch(hydroxyl_pattern):
        has_oxidation = True
    
    # Keto and aldehydes
    keto_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')
    if mol.HasSubstructMatch(keto_pattern):
        has_oxidation = True

    # Epoxide groups
    epoxide_pattern = Chem.MolFromSmarts('[OX2][CX4]1[OX2][CX4]1')
    if mol.HasSubstructMatch(epoxide_pattern):
        has_oxidation = True
    
    # Ester/Lactone groups
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    if mol.HasSubstructMatch(ester_pattern):
        has_oxidation = True

    # Carboxylic acid or related groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)O')
    if mol.HasSubstructMatch(carboxyl_pattern):
        has_oxidation = True

    if not has_oxidation:
        return False, "No typical oxidation patterns (hydroxy, keto, epoxy, ester) found"
    
    # Check for multiple unsaturations in carbon backbone
    cc_double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bond_matches = len(mol.GetSubstructMatches(cc_double_bond_pattern))
    if double_bond_matches < 1:
        return False, f"Expected at least 1 double bond, found {double_bond_matches}"

    # Check for other important structural motifs that may indicate analogues
    # Includes cyclic structures beyond just traditional rings
    num_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    if num_rings > 0:
        has_macrocycle = any(subs for subs in mol.GetSubstructMatches(Chem.MolFromSmarts('[R]')))
        if has_macrocycle:
            return True, "Contains both unsaturations and complex structural oxidation-responsive motifs characteristic of icosanoids"

    return True, "Contains an icosanoid-like backbone with oxidation and unsaturation"