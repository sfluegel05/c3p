"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins have a characteristic penam core (beta-lactam fused to thiazolidine ring)
    with specific substituents: two methyls at position 2, carboxylate at position 3,
    and carboxamido at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES with stereochemistry
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic penam core pattern - beta-lactam fused to thiazolidine
    penam_core = Chem.MolFromSmarts('[#7]1[#6][#6](=O)[#7][#6]2[#16][#6][#6]12')
    if penam_core is None:
        return None, "Invalid SMARTS pattern for penam core"
    if not mol.HasSubstructMatch(penam_core):
        return False, "No penam core structure found"

    # Check for the specific stereochemistry of the bicyclic system
    stereo_core = Chem.MolFromSmarts('[H][C@]12S[C@@]([C@@H](N1)C(=O)N2)(C)C')
    if stereo_core is None:
        return None, "Invalid SMARTS pattern for stereo core"
    if not mol.HasSubstructMatch(stereo_core, useChirality=True):
        return False, "Incorrect stereochemistry of penam core"

    # Check for carboxylate group (both acid and ester forms)
    carboxyl = Chem.MolFromSmarts('[CX3](=[OX1])[OX2,OX1-]')
    if carboxyl is None:
        return None, "Invalid SMARTS pattern for carboxyl"
    if not mol.HasSubstructMatch(carboxyl):
        return False, "Missing carboxylate group"

    # Check for two methyl groups at position 2
    dimethyl = Chem.MolFromSmarts('S[C](C)(C)')
    if dimethyl is None:
        return None, "Invalid SMARTS pattern for dimethyl"
    if not mol.HasSubstructMatch(dimethyl):
        return False, "Missing geminal dimethyl groups"

    # Check for carboxamido group at position 6
    carboxamido = Chem.MolFromSmarts('[NH][CH]1C(=O)N')
    if carboxamido is None:
        return None, "Invalid SMARTS pattern for carboxamido"
    if not mol.HasSubstructMatch(carboxamido):
        return False, "Missing carboxamido group"

    # Verify atom counts
    s_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#16]')))
    n_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#7]')))
    o_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#8]')))

    if s_count != 1:
        return False, f"Incorrect number of sulfur atoms (found {s_count}, expected 1)"
    if n_count < 2:
        return False, f"Too few nitrogen atoms (found {n_count}, expected at least 2)"
    if o_count < 3:
        return False, f"Too few oxygen atoms (found {o_count}, expected at least 3)"

    # Check ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Missing required bicyclic system"

    return True, "Structure contains penam core with correct substituents and stereochemistry"