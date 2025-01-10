"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: CHEBI:16480 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone with two fatty acid chains and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
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
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Check for phosphate group (P with 4 oxygens attached)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidic acids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for phosphatidic acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for phosphatidic acid"
    if o_count < 6:
        return False, "Too few oxygens for phosphatidic acid"

    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphate group"