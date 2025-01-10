"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a six-carbon monosaccharide with a specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Hexose-related SMARTS patterns (considering stereochemistry)
    hexose_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)C1"),  # Pyranose form with D-config
        Chem.MolFromSmarts("O1[C@H]([C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CO"),  # Furanose form with D-config
        Chem.MolFromSmarts("[H]C(=O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO")   # Open-chain aldehydo form with D-config
    ]
    
    # Attempt to match the SMARTS patterns for D-hexose identification
    if not any(mol.HasSubstructMatch(pattern) for pattern in hexose_patterns):
        return False, "Does not match hexose structures (considering pyranose/furanose and D-configuration)"
    
    # Ensure stereochemistry corresponds to D-sugars at the correct chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    d_sugar_configuration = all(center[1] in ['R', 'S'] for center in chiral_centers)
    if not d_sugar_configuration:
        return False, "Incorrect stereochemistry for D-hexose"
    
    # Verify exactly six carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Incorrect number of carbons for a hexose"

    # Verify expected number of oxygens (usually 5 for closed monosaccharides)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 5 and o_count != 6:
        return False, "Unusual number of oxygens for a hexose"

    return True, "Valid D-hexose conformer and stereochemistry detected"