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
    These compounds have a steroid core with a 3-beta hydroxyl group and a double bond between C5-C6.
    
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

    # Define a more flexible steroid core pattern that matches the tetracyclic system
    # This pattern allows for variations in ring connections and substitutions
    steroid_core = Chem.MolFromSmarts("[C,c]1[C,c][C,c][C,c]2[C,c]([C,c]1)[C,c][C,c]3[C,c][C,c][C,c]4[C,c][C,c][C,c]([C,c]4)[C,c]3[C,c]2")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 3-beta-hydroxyl group
    # [C@@H] specifies the beta stereochemistry at position 3
    # The pattern looks for the hydroxyl in the context of the A ring
    beta_hydroxyl_pattern = Chem.MolFromSmarts("[C,c]1[C,c][C,c][C@@H]([OH1])[C,c][C,c]1")
    
    if not mol.HasSubstructMatch(beta_hydroxyl_pattern):
        return False, "No 3-beta-hydroxyl group found"

    # Check for double bond between C5-C6
    # This pattern specifically looks for the double bond in ring B of the steroid
    delta_5_pattern = Chem.MolFromSmarts("[C,c]1[C,c][C,c]=C[C,c][C,c]1")
    
    if not mol.HasSubstructMatch(delta_5_pattern):
        return False, "No double bond between positions 5 and 6"

    # Additional structural checks
    
    # Count rings (steroids should have at least 4)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Count carbons (steroids typically have 17+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Count oxygens (should have at least one for the hydroxyl group)
    oxygen_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    if oxygen_count < 1:
        return False, "No oxygen atoms found"

    # Verify sp2 carbons (should have at least 2 for the double bond)
    sp2_carbons = len([atom for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP2])
    if sp2_carbons < 2:
        return False, "No carbon-carbon double bond found"

    return True, "Contains steroid core with 3-beta-hydroxyl group and double bond between C5-C6"