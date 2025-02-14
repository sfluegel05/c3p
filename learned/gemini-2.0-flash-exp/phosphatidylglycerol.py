"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol is a glycerol backbone with two fatty acid chains and a
    phosphoglycerol group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for two ester groups attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups on glycerol, need at least 2"

    # Look for the phosphatidyl group, which contains a phosphate diester connecting to a glycerol
    phospho_glycerol_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])[OX2][CH2X4][CHX4][CH2X4][OX2]")
    phospho_glycerol_matches = mol.GetSubstructMatches(phospho_glycerol_pattern)
    if len(phospho_glycerol_matches) < 1:
         return False, f"Found {len(phospho_glycerol_matches)} phosphoglycerol groups, expected at least 1"


    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"


    #Check rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - PGs are typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for phosphatidylglycerol"

    # Count carbons, oxygens and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if c_count < 15:
         return False, "Too few carbons for phosphatidylglycerol"
    if o_count < 8:
        return False, "Too few oxygens for phosphatidylglycerol"
    if p_count < 1:
        return False, "Must have at least one phosphorus"

    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphoglycerol group"