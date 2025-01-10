"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol has a glycerol backbone with a phosphatidyl group (phosphate + two fatty acids)
    and an additional glycerol moiety attached to the phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphate group (P with 4 oxygens, one of which is double-bonded)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for two ester groups (fatty acid chains attached via -O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for additional glycerol moiety attached to phosphate
    # This is a glycerol with a hydroxyl group and a phosphate group
    additional_glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2])[CH2X4][OX2][PX4]")
    if not mol.HasSubstructMatch(additional_glycerol_pattern):
        return False, "No additional glycerol moiety attached to phosphate"

    # Check for long carbon chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidylglycerols typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylglycerol"

    return True, "Contains glycerol backbone, phosphatidyl group, and additional glycerol moiety"