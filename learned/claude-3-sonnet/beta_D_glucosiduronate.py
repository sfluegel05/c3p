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
    Beta-D-glucosiduronate is the deprotonated form of beta-D-glucuronic acid
    conjugated to another molecule via a glycosidic bond.
    
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

    # First check: Must have carboxylate group
    carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate):
        return False, "No carboxylate group found"

    # Check for beta-D-glucuronate core with glycosidic linkage
    # [*] represents connection to rest of molecule
    # Specific stereochemistry pattern for beta-D-glucuronate:
    # - Beta configuration at anomeric carbon (C1)
    # - Hydroxyls at C2,C3,C4 with specific orientation
    # - Carboxylate at C5
    beta_glucosiduronate = Chem.MolFromSmarts(
        "[#6,#8,#7,#16;!$(*C(=O)[O-])][OX2][C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C(=O)[O-])O1"
    )
    
    if not mol.HasSubstructMatch(beta_glucosiduronate):
        return False, "No beta-D-glucosiduronate group found with correct stereochemistry and linkage"

    # Count atoms to ensure it's not just a simple sugar acid
    atom_count = mol.GetNumAtoms()
    if atom_count < 15:  # Beta-D-glucuronate alone would have ~12-13 atoms
        return False, "Molecule too small - appears to be standalone sugar acid"

    # Additional check for proper glycosidic linkage
    # The oxygen should be connected to another carbon/oxygen/nitrogen/sulfur
    glycosidic_link = Chem.MolFromSmarts("[#6,#8,#7,#16;!$(*C(=O)[O-])][OX2][C@H]1O[C@H][CH][CH]1")
    if not mol.HasSubstructMatch(glycosidic_link):
        return False, "No proper glycosidic linkage found"

    # Verify pyranose ring form
    pyranose = Chem.MolFromSmarts("O1[CH][CH][CH][CH][CH]1")
    if not mol.HasSubstructMatch(pyranose):
        return False, "Not in pyranose form"

    # Check for exactly 3 hydroxyls on carbons 2,3,4
    hydroxyls = Chem.MolFromSmarts("[CH]1([O][!C])[CH]([O][!C])[CH]([O][!C])[CH][CH](C(=O)[O-])O1")
    if not mol.HasSubstructMatch(hydroxyls):
        return False, "Incorrect hydroxyl pattern"

    return True, "Contains beta-D-glucosiduronate group with correct stereochemistry and glycosidic linkage"