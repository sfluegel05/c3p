"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:35741 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three substituent groups (acyl, alkyl, or alk-1-enyl) attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 3 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"

    # Check for substituent groups (acyl, alkyl, or alk-1-enyl) attached to the glycerol backbone
    # Acyl: -C(=O)-R
    # Alkyl: -C-R
    # Alk-1-enyl: -C=C-R
    substituent_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    if len(substituent_matches) < 3:
        return False, f"Missing substituent groups, got {len(substituent_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be substituent groups"

    # Check molecular weight - triradylglycerols typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for triradylglycerol"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for triradylglycerol"
    if o_count != 6:
        return False, "Must have exactly 6 oxygens (3 ester groups)"

    return True, "Contains glycerol backbone with 3 substituent groups attached via ester bonds"