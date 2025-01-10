"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    PIPs contain a myo-inositol ring with at least one phosphate group, connected via 
    phosphodiester linkage to a glycerol backbone with two fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PIP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring (cyclohexane with multiple OH groups)
    inositol_pattern = Chem.MolFromSmarts("[CH1]1[CH1][CH1][CH1][CH1][CH1]1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for phosphate groups (-OP(=O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # Look for glycerol backbone (C-C-C with O attachments)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Count phosphorus atoms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 1:
        return False, "No phosphorus atoms found"

    # Count oxygen atoms (PIPs typically have many oxygens from phosphates and OH groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 10:
        return False, "Too few oxygen atoms for a PIP"

    # Verify presence of fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Check molecular weight (PIPs are large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for PIP"

    # Look for phosphodiester linkage between glycerol and inositol
    phosphodiester_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage found"

    return True, "Contains inositol ring with phosphate(s), glycerol backbone with fatty acids, and phosphodiester linkage"