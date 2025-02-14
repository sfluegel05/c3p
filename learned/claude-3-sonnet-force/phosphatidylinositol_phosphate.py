"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: CHEBI:28085 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 1 or 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X3,CH2X4][CHX4][CH2X3,CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for at least 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Look for inositol ring with 1 or more phosphate groups
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)OP(=O)([O-,O])[O-])OP(=O)([O-,O])[O-])OP(=O)([O-,O])[O-])OP(=O)([O-,O])[O-])OP(=O)([O-,O])[O-]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring with phosphate groups found"

    # Look for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidylinositol phosphates typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for phosphatidylinositol phosphate"

    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphorylated inositol head group"