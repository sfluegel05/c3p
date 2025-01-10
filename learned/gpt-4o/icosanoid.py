"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are characterized by essential fatty acid oxidation products often having C20 backbones with unsaturation and oxidation sites.

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
    
    # Look for carbon backbone typical of EFAs, extended range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 24:
        return False, f"Expected 18-24 carbons, found {c_count}"
    
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
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    if mol.HasSubstructMatch(epoxide_pattern):
        has_oxidation = True
    
    # Ester/Lactone groups
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    if mol.HasSubstructMatch(ester_pattern):
        has_oxidation = True

    if not has_oxidation:
        return False, "No typical oxidation patterns (hydroxy, keto, epoxy, ester) found"
    
    # Check for multiple unsaturations in carbon backbone
    cc_double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bond_matches = len(mol.GetSubstructMatches(cc_double_bond_pattern))
    if double_bond_matches < 2:
        return False, f"Expected at least 2 double bonds, found {double_bond_matches}"
    
    # Check for other important structural motifs, such as:
    # - Cyclic structures or specific ring types present in analogues
    # - Core structures that define sub-classes

    num_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    if (num_rings > 0):
        additional_pattern = Chem.MolFromSmarts('[R]')
        if mol.HasSubstructMatch(additional_pattern):
            return True, "Additional structural patterns matched"

    # If all patterns match
    return True, "Contains an icosanoid-like C18-24 backbone with oxidation and multiple unsaturations"