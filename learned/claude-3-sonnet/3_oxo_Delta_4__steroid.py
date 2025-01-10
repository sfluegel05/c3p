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

    # Check for steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for 3-oxo group and Delta-4 double bond pattern
    # This SMARTS pattern looks for:
    # - A ketone (C=O) at position 3
    # - A double bond between carbons 4 and 5
    # - The correct connectivity in the A-ring of the steroid
    oxo_delta4_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6](=O)~[#6]=[#6]~[#6]~1")
    
    if not mol.HasSubstructMatch(oxo_delta4_pattern):
        return False, "No 3-oxo-Delta(4) pattern found"

    # Additional check to ensure the ketone and double bond are conjugated
    conjugated_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6]=[#6]")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Ketone and double bond not conjugated"

    # Count carbons to ensure reasonable size for a steroid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19:  # Most steroids have at least 19 carbons
        return False, "Too few carbons for a steroid structure"

    # Check ring count (steroids typically have 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    return True, "Contains steroid core with 3-oxo group conjugated to Delta-4 double bond"