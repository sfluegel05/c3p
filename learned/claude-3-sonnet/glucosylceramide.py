"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Glucosylceramides are cerebrosides with a glucose head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for beta-D-glucose moiety
    # Pattern matches beta-D-glucose with correct stereochemistry at all centers
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([OH0,OH1])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check for ceramide core structure with more flexible matching
    # Allow for variations in sphingoid base and fatty acid chain
    ceramide_core = Chem.MolFromSmarts("[CH2,CH][CH2,CH][CH2,CH][#6]~[#6]~[#6][CH]([OH1])[CH]([NH1][C](=O))[CH2]O")
    if not mol.HasSubstructMatch(ceramide_core):
        return False, "No ceramide core structure found"

    # Check for glycosidic linkage between glucose and ceramide
    # More flexible pattern that allows for variations in linkage geometry
    glycosidic_link = Chem.MolFromSmarts("[OH0][CH2][CH]([NH1])[CH]([OH1])")
    if not mol.HasSubstructMatch(glycosidic_link):
        return False, "No glycosidic linkage found"

    # Check for long chain fatty acid
    fatty_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatches(fatty_chain):
        return False, "No long fatty acid chain found"

    # Count atoms for basic composition check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Basic composition requirements
    if c_count < 20:  # Minimum carbons for smallest glucosylceramide
        return False, "Insufficient carbon atoms for glucosylceramide"
    if o_count < 7:  # Minimum oxygens needed
        return False, "Insufficient oxygen atoms"
    if n_count != 1:  # Must have exactly one nitrogen in amide bond
        return False, "Must have exactly one nitrogen atom"

    # Additional structural checks
    # Verify presence of amide group
    amide_pattern = Chem.MolFromSmarts("[NH1][C](=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Verify hydroxyl groups on sphingoid base
    sphingoid_hydroxyls = Chem.MolFromSmarts("[CH]([OH1])[CH]([NH1])")
    if not mol.HasSubstructMatch(sphingoid_hydroxyls):
        return False, "Missing characteristic hydroxyl groups"

    return True, "Contains beta-D-glucose linked to ceramide with appropriate structural features"