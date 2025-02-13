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

    # Check for basic steroid core (four fused rings)
    # More flexible pattern that allows for variations
    steroid_core = Chem.MolFromSmarts("C1C[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # Pattern specifically looks for C=O at position 3 of the A ring
    oxo_3_pattern = Chem.MolFromSmarts("[CH2][CH2]C(=O)[CH2][C@@H]")
    if not mol.HasSubstructMatch(oxo_3_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-alpha configuration
    # In 5-alpha steroids, rings A/B are trans-fused
    # The hydrogen at C5 is in alpha (below plane) orientation
    alpha_5_pattern = Chem.MolFromSmarts("[C]1[CH2]C(=O)[CH2][C@@H]([CH2])[CH2]")
    if not mol.HasSubstructMatch(alpha_5_pattern):
        return False, "No 5-alpha configuration found"

    # Verify A/B ring fusion stereochemistry
    ab_fusion_pattern = Chem.MolFromSmarts("[C]1[CH2]C(=O)[CH2][C@@H]2[CH2][CH2]")
    if not mol.HasSubstructMatch(ab_fusion_pattern):
        return False, "Incorrect A/B ring fusion stereochemistry"

    # Additional validation checks
    
    # Count carbons (steroids typically have 19+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19:
        return False, "Too few carbons for a steroid structure"

    # Check for reasonable molecular weight
    mol_wt = sum([atom.GetMass() for atom in mol.GetAtoms()])
    if mol_wt < 250 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for steroids"

    # Check ring count (steroids should have at least 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    return True, "Molecule contains 3-oxo-5alpha-steroid structure with correct stereochemistry"