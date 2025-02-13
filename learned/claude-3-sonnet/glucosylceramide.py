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

    # Check for beta-D-glucose moiety with specific stereochemistry
    # This pattern enforces the exact stereochemistry of beta-D-glucose
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([OH0,OH1])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Verify it's not galactose (different C4 stereochemistry)
    galactose_pattern = Chem.MolFromSmarts("[C@@H]1([OH0,OH1])O[C@H](CO)[C@@H](O)[C@@H](O)[C@H]1O")
    if mol.HasSubstructMatch(galactose_pattern):
        return False, "Contains galactose instead of glucose"

    # Check for ceramide core with flexible matching
    # Allow for various sphingoid base modifications
    ceramide_core = Chem.MolFromSmarts("[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2][CH]([OH1])[CH]([NH1][C](=O))[CH2]O")
    if not mol.HasSubstructMatch(ceramide_core):
        return False, "No ceramide core structure found"

    # Verify glycosidic linkage to ceramide
    glycosidic_link = Chem.MolFromSmarts("[C@@H]1([OH0]C[C@H]([NH1])[C@H](O))O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycosidic_link):
        return False, "Incorrect glycosidic linkage"

    # Check for fatty acid chain
    fatty_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No long fatty acid chain found"

    # Count atoms for composition check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    # Verify composition
    if c_count < 24:  # Minimum carbons for smallest glucosylceramide
        return False, "Insufficient carbon atoms for glucosylceramide"
    if o_count < 8 or o_count > 9:  # Glucose (6) + ceramide (2-3)
        return False, "Incorrect number of oxygen atoms"
    if n_count != 1:  # Must have exactly one nitrogen in amide bond
        return False, "Must have exactly one nitrogen atom"
    if p_count > 0:  # Exclude phosphorylated derivatives
        return False, "Contains phosphate group"

    return True, "Contains beta-D-glucose linked to ceramide with appropriate structural features"