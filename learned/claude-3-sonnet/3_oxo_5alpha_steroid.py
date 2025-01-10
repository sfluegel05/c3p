"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid has:
    - A steroid core (four fused rings)
    - A ketone group at position 3
    - Alpha (below plane) configuration at position 5
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - more flexible version
    # Matches the four-ring system without being too specific about substituents
    steroid_core = Chem.MolFromSmarts("[C]1[C][C]2[C][C][C]3[C][C][C]4[C][C][C][C]4[C]3[C][C]2[C]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 3-oxo group
    # Pattern looks for ketone at position 3, allowing for various substituents
    oxo_3_pattern = Chem.MolFromSmarts("[CH2][CH2]C(=O)[CH2][CH]")
    if not mol.HasSubstructMatches(oxo_3_pattern):
        return False, "No ketone group at position 3"

    # Check 5-alpha stereochemistry
    # This pattern specifically looks for the trans fusion of rings A/B
    # with the hydrogen at C5 in alpha (below plane) orientation
    alpha_5_pattern = Chem.MolFromSmarts("[C]1[CH2]C(=O)[CH2][C@@H]([CH2])[CH2]")
    if not mol.HasSubstructMatch(alpha_5_pattern):
        return False, "No 5-alpha configuration found"

    # Additional validation of ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (most steroids have at least 19 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:  # Being more lenient with minimum carbon count
        return False, "Too few carbons for steroid structure"

    # Check for reasonable molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:  # More lenient range
        return False, "Molecular weight outside typical range for steroids"

    # Verify basic connectivity of steroid core
    # Look for the characteristic four-ring system with correct fusion points
    steroid_connectivity = Chem.MolFromSmarts("C1CC2CCC3C4CCCC4CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_connectivity):
        return False, "Incorrect steroid ring connectivity"

    # Additional check for A/B ring trans fusion characteristic of 5-alpha steroids
    ab_fusion = Chem.MolFromSmarts("[C]1[CH2]C(=O)[CH2][C@@H]2[CH2][CH2]")
    if not mol.HasSubstructMatch(ab_fusion):
        return False, "Incorrect A/B ring fusion stereochemistry"

    return True, "Contains 3-oxo-5alpha-steroid structure with correct stereochemistry"