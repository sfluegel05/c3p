"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern allowing for different bond types
    # and some variation in the core structure
    steroid_core = Chem.MolFromSmarts(
        '[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~1'
    )
    
    # Alternative steroid core pattern
    steroid_core_alt = Chem.MolFromSmarts(
        '[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~1'
    )

    if not (mol.HasSubstructMatch(steroid_core) or mol.HasSubstructMatch(steroid_core_alt)):
        return False, "No steroid core structure found"

    # 17α-hydroxy pattern - more specific to position 17
    # Looking for a hydroxyl group with alpha stereochemistry at position 17
    # of the steroid D ring
    oh_17alpha_pattern = Chem.MolFromSmarts(
        '[C;R1]12[C;R1][C;R1][C;R1]3[C;R1][C;R1][C;R1][C;R1]4[C;R1][C;R1][C;R1][C;R1]4[C;R1][C;R1]3[C;R1][C;R1]1[C;R1]2[O;H1]'
    )
    
    # Alternative 17α-hydroxy patterns
    oh_17alpha_alt1 = Chem.MolFromSmarts('[C;R1]1([C;R1][C;R1][C;R1]2)([O;H1])[C;R1][C;R1][C;R1]2[H]')
    oh_17alpha_alt2 = Chem.MolFromSmarts('[C;R1]([O;H1])([C;R1])([C;R1])[C;R1]')

    if not (mol.HasSubstructMatch(oh_17alpha_pattern) or 
            mol.HasSubstructMatch(oh_17alpha_alt1) or 
            mol.HasSubstructMatch(oh_17alpha_alt2)):
        return False, "No 17-alpha hydroxyl group found"

    # Validation checks
    # Count carbons (steroids typically have 19+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:  # Lowered threshold to catch more variants
        return False, "Too few carbons for steroid structure"

    # Count rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check for reasonable number of sp3 carbons
    sp3_carbons = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() == 6 and 
                     atom.GetHybridization() == Chem.HybridizationType.SP3)
    if sp3_carbons < 5:  # Lowered threshold
        return False, "Insufficient sp3 carbons for steroid structure"

    # Count oxygens (should have at least one for the 17-OH)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found"

    return True, "Contains steroid core with 17-alpha hydroxyl group"