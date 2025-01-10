"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is a cytidine diphosphate (CDP) attached to a glycerol backbone with two acyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cytidine diphosphate (CDP) group
    cdp_pattern = Chem.MolFromSmarts("[NX3]C1=CC(=O)[NX3]C(=N1)[CX4][CX4][OX2][PX4](=[OX1])([OX2][PX4](=[OX1])[OX1])[OX1]")
    if cdp_pattern is None or not mol.HasSubstructMatch(cdp_pattern):
        return False, "No cytidine diphosphate (CDP) group found"

    # Check for glycerol backbone with two acyl groups
    glycerol_pattern = Chem.MolFromSmarts("[CX4][CX4]([OX2][CX4](=[OX1])[CX4])[CX4]([OX2][CX4](=[OX1])[CX4])[CX4][OX2][PX4](=[OX1])([OX2][PX4](=[OX1])[OX1])[OX1]")
    if glycerol_pattern is None or not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two acyl groups found"

    # Check for two acyl groups (ester bonds)
    ester_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4](=[OX1])[CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl groups, need at least 2"

    # Check for long carbon chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - CDP-diacylglycerols typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for CDP-diacylglycerol"

    return True, "Contains cytidine diphosphate (CDP) group, glycerol backbone, and two acyl groups"