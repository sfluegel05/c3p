"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three radyl groups attached.
    Radyl groups can be acyl, alkyl or alk-1-enyl.

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

    # 1. Glycerol backbone check: same as for triglycerides.
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # 2. Check for three radyl groups - we combine the SMARTS of possible linkages
    radyl_pattern = Chem.MolFromSmarts("[OX2]-[CX3,CX4]=[OX1];!@[CH2X4][CHX4][CH2X4], [OX2]-[CX4];!@[CH2X4][CHX4][CH2X4],[OX2]-[CX3]=[CX3];!@[CH2X4][CHX4][CH2X4]")
    radyl_matches = mol.GetSubstructMatches(radyl_pattern)

    if len(radyl_matches) != 3:
       return False, f"Found {len(radyl_matches)} radyl groups, need exactly 3"

    # 3. Check for long carbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbons ({c_count}) in chains"

    # 4. Check number of oxygens: should be 6
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 6:
        return False, "Must have exactly 6 oxygens"
    
    # 5. Check presence of other elements, apart from C, H and O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1,6,8]:
            return False, "Molecule contains elements other than C, H and O"

    return True, "Contains glycerol backbone with three radyl groups (acyl/alkyl/alk-1-enyl)"