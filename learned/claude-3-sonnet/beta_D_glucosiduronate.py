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

    # Basic pyranose ring with carboxylate and any glycosidic linkage
    # Note: [#6,#8,#7,#16] represents C,O,N,S atoms that could be part of the rest of molecule
    basic_structure = Chem.MolFromSmarts("""
        [#6,#8,#7,#16]O[CH]1[CH]([OH])[CH]([OH])[CH]([OH])[CH](C(=O)[O-])O1
    """)
    
    if not mol.HasSubstructMatch(basic_structure):
        return False, "No glucuronate core structure found"

    # Check stereochemistry at key carbons
    # Beta configuration at anomeric carbon (C1)
    beta_config = Chem.MolFromSmarts("""
        [#6,#8,#7,#16]O[C@H]1O[C@H](C(=O)[O-])[C@@H]([OH])[C@H]([OH])[C@@H]1[OH]
    """)
    
    if not mol.HasSubstructMatch(beta_config):
        return False, "Incorrect stereochemistry - must be beta configuration"

    # Count atoms to ensure it's not just a simple sugar acid
    atom_count = mol.GetNumAtoms()
    if atom_count < 15:  # Beta-D-glucuronate alone would have ~12-13 atoms
        return False, "Molecule too small - appears to be standalone sugar acid"

    # Additional check for proper glycosidic linkage
    # Look for oxygen connected to the anomeric carbon and another heavy atom
    glycosidic = Chem.MolFromSmarts("""
        [#6,#8,#7,#16;!$(*C(=O)[O-])]O[C@H]1O[C@H](C(=O)[O-])
    """)
    
    if not mol.HasSubstructMatch(glycosidic):
        return False, "No proper glycosidic linkage found"

    # Verify hydroxyl groups at C2, C3, C4 positions
    hydroxyls = Chem.MolFromSmarts("""
        [OH][CH]-[CH]([OH])-[CH]([OH])-[CH]-[CH](C(=O)[O-])O
    """)
    
    if not mol.HasSubstructMatch(hydroxyls):
        return False, "Missing required hydroxyl groups"

    return True, "Contains beta-D-glucosiduronate group with correct stereochemistry and glycosidic linkage"