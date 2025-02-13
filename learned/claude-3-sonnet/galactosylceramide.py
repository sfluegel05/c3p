"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide has a galactose sugar connected to a ceramide 
    (sphingosine + fatty acid) via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pyranose ring pattern (more general than previous version)
    # Matches 6-membered sugar ring with oxygen
    sugar_ring = Chem.MolFromSmarts("C1OCCCC1")
    if not mol.HasSubstructMatch(sugar_ring):
        return False, "No sugar ring found"

    # Look for galactose hydroxyl pattern - at least 3 OH groups on ring carbons
    # and one CH2OH group
    galactose_pattern = Chem.MolFromSmarts("[CH2][OH]")
    oh_pattern = Chem.MolFromSmarts("[CH]([OH])")
    
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "Missing CH2OH group characteristic of galactose"
    
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, "Insufficient hydroxyl groups for galactose"

    # Look for amide bond (-NH-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")
    if not mol.HasSubstructMatches(amide_pattern):
        return False, "No amide bond found"

    # Look for long carbon chains (sphingosine base and fatty acid)
    # More flexible pattern to match various chain lengths
    carbon_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2]")
    chain_matches = len(mol.GetSubstructMatches(carbon_chain))
    if chain_matches < 2:
        return False, "Missing required long carbon chains"

    # Look for glycosidic linkage
    glycosidic_pattern = Chem.MolFromSmarts("[CH2]O[CH]1O[CH][CH][CH][CH]C1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Check for characteristic OH groups of sphingosine/ceramide
    sphingosine_pattern = Chem.MolFromSmarts("[CH]([OH])[CH]([NH])")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Missing characteristic OH groups of ceramide"

    # Count atoms to ensure molecule is large enough
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:  # Adjusted minimum carbon count
        return False, "Too few carbons for galactosylceramide"
    if o_count < 6:  # Minimum oxygen count for basic structure
        return False, "Too few oxygens for galactosylceramide"

    # Check for optional sulfate modification
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[OH]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    base_reason = "Contains galactose connected to ceramide via glycosidic bond"
    if has_sulfate:
        return True, base_reason + " with sulfate modification"
    return True, base_reason