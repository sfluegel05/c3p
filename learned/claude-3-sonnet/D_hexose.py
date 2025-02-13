"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16234 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    D-hexose is a hexose (6-carbon monosaccharide) with D-configuration at C5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular formula (allowing for isotope variations like 13C)
    base_atoms = {'C': 6, 'H': 12, 'O': 6}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in base_atoms:
            return False, f"Contains unexpected element: {symbol}"
        if symbol == 'C' and atom.GetIsotope() not in [0, 13]:
            return False, "Contains unexpected carbon isotope"

    # Define SMARTS patterns for different hexose forms
    # Pyranose patterns with D-configuration at C5
    d_pyranose_patterns = [
        # alpha-D-pyranose pattern
        "[C@@H]1(CO)O[CH](O)[CH]([CH]([CH]([CH]1O)O)O)O",
        # beta-D-pyranose pattern
        "[C@H]1(CO)O[CH](O)[CH]([CH]([CH]([CH]1O)O)O)O",
        # Open-chain D-aldose pattern with explicit D-configuration at C5
        "O=C[CH](O)[CH](O)[CH](O)[C@H](O)CO",
        # D-furanose patterns
        "O1[CH]([C@H](O)CO)[CH](O)[CH](O)[CH]1O"
    ]

    # Convert patterns to RDKit molecules
    pattern_mols = [Chem.MolFromSmarts(p) for p in d_pyranose_patterns]
    
    # Check if molecule matches any D-hexose pattern
    matches_any = False
    for pattern in pattern_mols:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matches_any = True
            break
            
    if not matches_any:
        return False, "Does not match D-hexose pattern"

    # Additional checks for correct structure
    # Count chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol)
    if len(chiral_centers) < 3:
        return False, f"Too few chiral centers for hexose: {len(chiral_centers)}"

    # Count hydroxy groups
    oh_pattern = Chem.MolFromSmarts("[OH]")
    oh_count = len(mol.GetSubstructMatches(oh_pattern))
    if oh_count < 4:
        return False, "Too few hydroxyl groups for hexose"

    # Check for carbonyl or hemiacetal carbon
    carbonyl_pattern = Chem.MolFromSmarts("[$([CH]=O),$([CH]1O[CH](O)[CH]([CH]([CH]([CH]1O)O)O)O)]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing required carbonyl or hemiacetal group"

    # Verify carbon chain connectivity
    chain_pattern = Chem.MolFromSmarts("CC(O)C(O)C(O)C(O)CO")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Incorrect carbon chain connectivity"

    return True, "Matches D-hexose pattern with correct stereochemistry"