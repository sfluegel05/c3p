"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:35353 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    These steroids have a ketone at position 3 conjugated with a C=C double bond at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings (A,B,C,D)
    # More flexible pattern that focuses on the connectivity of the rings
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~[#6]~[#6]~1"
    )
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Pattern for 3-oxo group and Delta-4 double bond in ring A
    # This pattern specifically looks for:
    # - Ketone (C=O) at position 3
    # - Double bond between carbons 4 and 5
    # - Proper connectivity to the B ring
    oxo_delta4_pattern = Chem.MolFromSmarts(
        "[#6]~1~[#6]~C(=O)~C=C~[#6]~2~[#6]~[#6]~[#6]~1~[#6]~2"
    )
    
    if not mol.HasSubstructMatch(oxo_delta4_pattern):
        return False, "No 3-oxo-Delta(4) pattern found in ring A"

    # Additional check for conjugation between ketone and double bond
    conjugated_pattern = Chem.MolFromSmarts("C(=O)C=C")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Ketone not conjugated with double bond"

    # Validate basic steroid characteristics
    # Count carbons (steroids typically have 19-35 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (19 <= carbon_count <= 40):  # Increased upper limit to accommodate larger derivatives
        return False, f"Invalid carbon count ({carbon_count}) for steroid structure"

    # Count rings (steroids must have at least 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check that the molecule is not too small or too large
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if not (250 <= mol_weight <= 1000):  # Typical range for steroids and their derivatives
        return False, f"Molecular weight {mol_weight} outside typical steroid range"

    return True, "Contains steroid core with 3-oxo group conjugated to Delta-4 double bond"