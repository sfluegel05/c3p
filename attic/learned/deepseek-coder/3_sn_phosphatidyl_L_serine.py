"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:61811 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine is a glycerophosphoserine compound with acyl substituents at the 1- and 2-hydroxy positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a glycerol backbone with a phosphoserine group
    # The pattern should match: [glycerol]-[phosphate]-[serine]
    phosphoserine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]([OX2][PX4](=[OX1])([OX2][CH2X4][CHX4]([NX3])[CX3](=[OX1])[OX2H]))")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No glycerol backbone with phosphoserine group found"

    # Check for two ester groups (acyl chains) attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CX4][OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for the presence of serine moiety
    serine_pattern = Chem.MolFromSmarts("[NX3][CHX4][CX3](=[OX1])[OX2H]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine moiety found"

    # Check for long carbon chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - 3-sn-phosphatidyl-L-serine typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for 3-sn-phosphatidyl-L-serine"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for 3-sn-phosphatidyl-L-serine"
    if o_count < 6:
        return False, "Too few oxygens for 3-sn-phosphatidyl-L-serine"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (from serine)"

    return True, "Contains glycerol backbone with phosphoserine group and two acyl chains at the 1- and 2-positions"