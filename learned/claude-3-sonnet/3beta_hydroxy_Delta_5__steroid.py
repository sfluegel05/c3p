"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (four connected rings)
    # More flexible pattern that matches the four-ring system
    steroid_core = Chem.MolFromSmarts("C1C2CCC3CCCC4CCCC(C4)C3C2CC1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 3-beta-hydroxyl group with correct stereochemistry
    # This pattern specifically looks for the 3-beta-OH configuration in ring A
    # [C@@H] ensures the beta stereochemistry
    beta_hydroxyl_3 = Chem.MolFromSmarts("[CH2][C@@H](O)[CH2]")
    if not mol.HasSubstructMatch(beta_hydroxyl_3):
        return False, "No 3-beta-hydroxyl group found"

    # Check for double bond between C5-C6 (in ring B)
    # Pattern looks for the specific double bond in the context of the ring system
    double_bond_5_6 = Chem.MolFromSmarts("C1CC=CC2C1")
    if not mol.HasSubstructMatch(double_bond_5_6):
        return False, "No double bond between positions 5 and 6"

    # Additional check for the complete 3beta-hydroxy-Delta5 pattern
    # This combines the hydroxyl group and double bond in the correct positions
    complete_pattern = Chem.MolFromSmarts("[CH2][C@@H](O)C[CH2]C=C")
    if not mol.HasSubstructMatch(complete_pattern):
        return False, "Pattern of 3beta-OH and Delta-5 double bond not found in correct arrangement"

    # Verify carbon count (steroids typically have 17+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Check for reasonable number of rings (steroids should have 4 or more rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    return True, "Contains steroid core with 3-beta-hydroxyl group and double bond between C5-C6"