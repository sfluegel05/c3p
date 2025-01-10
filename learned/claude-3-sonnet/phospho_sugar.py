"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:37600 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide with a phosphate group esterified to a hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight bounds (150-500 Da for phospho sugars)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 150:
        return False, "Molecule too small to be a phospho sugar"
    if mol_weight > 500:
        return False, "Molecule too large to be a simple phospho sugar"

    # Check for excessive complexity
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:
        return False, "Too many rings for a simple phospho sugar"

    # Look for phosphate group
    phosphate_patterns = [
        "[OX2][P](=[O])([O,OH])[O,OH]",  # Any phosphate group
    ]
    
    has_phosphate = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_phosphate = True
            break
            
    if not has_phosphate:
        return False, "No phosphate group found"

    # Sugar patterns (more flexible now)
    sugar_patterns = [
        # Pyranose (6-membered ring) with optional deoxy
        "[CH2,O][CH1]1[CH1][CH1][CH1][CH1]O1",
        # Furanose (5-membered ring) with optional deoxy
        "[CH2,O][CH1]1[CH1][CH1][CH1]O1",
        # Open chain aldose/ketose
        "[CH2][CH1](-[O,H])[CH1](-[O,H])[CH1](-[O,H])[CH1]=O",
        "[CH2](-[O,H])[C](-[O,H])[CH1](-[O,H])[CH1](-[O,H])[CH2][O,H]",
        # Deoxy sugar patterns
        "[CH2]1[CH2][CH1]([OH])[CH1]([OH])O1",
        "[CH2]1[CH1]([OH])[CH1]([OH])[CH2]O1"
    ]

    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar = True
            break

    if not has_sugar:
        return False, "No sugar moiety found"

    # Exclude nucleotides and nucleosides
    nucleobase_patterns = [
        "c1ncnc2[nH]cnc12",  # Purine
        "c1c[nH]c(=O)[nH]c1=O",  # Pyrimidine
    ]
    
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Molecule appears to be a nucleotide/nucleoside"

    # Verify phosphate is connected to sugar
    sugar_phosphate_pattern = Chem.MolFromSmarts("[CH2,CH1]-[OX2]-[P](=[O])([O,OH])[O,OH]")
    if not mol.HasSubstructMatch(sugar_phosphate_pattern):
        return False, "Phosphate not directly connected to sugar carbon"

    # Count carbons to ensure it's a monosaccharide
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count > 8:
        return False, "Too many carbons for a simple phospho sugar"

    return True, "Contains monosaccharide with phosphate ester"