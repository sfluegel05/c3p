"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucosiduronate group.
    Beta-D-glucosiduronate is the deprotonated form of beta-D-glucuronic acid,
    with specific stereochemistry at all centers.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-glucuronate specifically
    # This pattern enforces:
    # - Beta configuration at C1 (DOWN bond)
    # - Specific stereochemistry at C2,C3,C4 (all DOWN except C3 UP)
    # - Carboxylate at C5
    # - Proper pyranose ring structure
    beta_glucuronate = Chem.MolFromSmarts(
        "[OX2][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C([O-])=O)O1"
    )
    
    if not mol.HasSubstructMatch(beta_glucuronate):
        return False, "No beta-D-glucuronate group found with correct stereochemistry"

    # Additional check for pyranose ring to exclude furanose forms
    pyranose_ring = Chem.MolFromSmarts("O1[CH][CH][CH][CH][CH]1")
    if not mol.HasSubstructMatch(pyranose_ring):
        return False, "No pyranose ring structure found"

    # Verify carboxylate is properly positioned
    carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate):
        return False, "No carboxylate group found"

    # Check for exactly 3 hydroxyl groups on the sugar ring
    ring_hydroxyls = Chem.MolFromSmarts("[CH]1([OH])[CH]([OH])[CH]([OH])[CH][CH](C(=O)[O-])O1")
    if not mol.HasSubstructMatch(ring_hydroxyls):
        return False, "Incorrect hydroxyl pattern on sugar ring"

    # Additional check for beta configuration
    # The O-glycosidic bond should be equatorial (beta)
    beta_config = Chem.MolFromSmarts("[OX2][C@H]1O[C@H][CH][CH][CH]1")
    if not mol.HasSubstructMatch(beta_config):
        return False, "Not in beta configuration"

    # Verify it's a glucuronate (not galacturonate or other sugar acids)
    # by checking specific stereochemistry pattern
    glucose_stereo = Chem.MolFromSmarts("[OX2][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C(=O)[O-])O1")
    if not mol.HasSubstructMatch(glucose_stereo):
        return False, "Not a glucuronate stereoisomer"

    return True, "Contains beta-D-glucuronate group with correct stereochemistry and carboxylate"