"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone - any oxolane having an oxo- substituent at any position 
on the tetrahydrofuran ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone contains a tetrahydrofuran ring with a ketone group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 5-membered ring containing oxygen (tetrahydrofuran/oxolane)
    # [O,o] means oxygen in or not in ring
    # ~1~1 means connected in a ring of size 5
    thf_pattern = Chem.MolFromSmarts("[O,o]~1~[#6]~[#6]~[#6]~[#6]~1")
    if not mol.HasSubstructMatch(thf_pattern):
        return False, "No tetrahydrofuran ring found"

    # Look for ketone group (C=O) 
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"
    
    # Two main types of tetrahydrofuranones:
    
    # 1. Lactone pattern (ketone is part of the ring)
    lactone_pattern = Chem.MolFromSmarts("O1CC[CH2,CH1,C]C(=O)1")
    
    # 2. Tetrahydrofuran with ketone substituent
    # This pattern looks for a tetrahydrofuran ring where one carbon has a C=O substituent
    thf_ketone_pattern = Chem.MolFromSmarts("[O;R1]1[CH2,CH1,C][CH2,CH1,C][CH2,CH1,C]([CH2,CH1,C]1)C(=O)")
    
    if not (mol.HasSubstructMatch(lactone_pattern) or mol.HasSubstructMatch(thf_ketone_pattern)):
        return False, "No ketone group attached to tetrahydrofuran ring"

    # Additional check to ensure no aromatic rings in the core structure
    # (tetrahydrofuran should be saturated)
    aromatic_thf = Chem.MolFromSmarts("[o]1cccc1")
    if mol.HasSubstructMatch(aromatic_thf):
        return False, "Core ring structure must be saturated (tetrahydro-)"

    return True, "Contains tetrahydrofuran ring with ketone group"