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

    # Basic steroid core pattern - four fused rings
    # More flexible pattern that captures the essential connectivity
    steroid_core = Chem.MolFromSmarts("C1C[C@H]2[C@H]3CC[C@H]4CCCC4[C@H]3CC[C@]12C")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 3-oxo group
    # Look for ketone at position 3 in ring A
    oxo_3_pattern = Chem.MolFromSmarts("[CH2][CH2]C(=O)[CH2][C@@H]")
    if not mol.HasSubstructMatch(oxo_3_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-alpha configuration
    # This looks for the characteristic trans fusion between rings A and B
    # with the hydrogen at C5 in alpha (below plane) orientation
    alpha_5_patterns = [
        Chem.MolFromSmarts("[CH2][CH2]C(=O)[CH2][C@@H]([CH2])[CH2]"),  # Basic pattern
        Chem.MolFromSmarts("[CH2][CH2]C(=O)[CH2][C@@H]([CH2])C"),      # With substitution
    ]
    
    has_5alpha = False
    for pattern in alpha_5_patterns:
        if mol.HasSubstructMatch(pattern):
            has_5alpha = True
            break
    
    if not has_5alpha:
        return False, "No 5-alpha configuration found"

    # Verify the presence of four fused rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Basic size check
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:  # Most steroids have at least 20 atoms
        return False, "Molecule too small for steroid structure"

    return True, "Contains 3-oxo-5alpha-steroid structure with correct stereochemistry"