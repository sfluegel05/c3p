"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate has a glycerol backbone, a phosphate group,
    a myo-inositol ring, and fatty acid chains.

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

    # 1. Check for glycerol backbone with O attachments (O-C-C-C-O)
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CH2X4][CHX4][CH2X4][OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
         return False, "No glycerol backbone found"


    # 2. Check for phosphate group connected to glycerol and inositol. Look for -C-O-P(=O)(O)-O-Inositol
    phosphate_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[PX4](=[OX1])(-[OX2])-[OX2]-[C]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate connecting glycerol and inositol"

    # 3. Check for myo-inositol ring (six-membered ring with 6 hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("C1C(C(C(C(C1O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for different myo-inositol phosphate substitutions
    # PIP can have 1, 2 or 3 phosphates in positions 3, 4 and 5 on the inositol ring
    has_inositol_phosphate = False
    inositol_phosphate_patterns = [
        Chem.MolFromSmarts("C1[C@H]([C@H]([C@H]([C@@H]([C@H]1O)O)O)OP(=O)(O)O)"), # 4-phosphate
        Chem.MolFromSmarts("C1[C@H]([C@H]([C@H]([C@@H]([C@H]1O)OP(=O)(O)O)O)O)"),  # 5-phosphate
        Chem.MolFromSmarts("C1[C@H]([C@@H]([C@H]([C@@H]([C@H]1O)OP(=O)(O)O)OP(=O)(O)O)O)"), # 4,5-bisphosphate
         Chem.MolFromSmarts("C1[C@H]([C@H]([C@@H]([C@@H]([C@H]1O)O)OP(=O)(O)O)OP(=O)(O)O)"), # 3,4-bisphosphate
        Chem.MolFromSmarts("C1[C@H]([C@@H]([C@H]([C@@H]([C@H]1O)OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)"), # 3,4,5 trisphosphate
        Chem.MolFromSmarts("C1[C@H]([C@@H]([C@@H]([C@H]([C@H]1O)O)OP(=O)(O)O)O)"), # 3-phosphate
      
    ]
    for pattern in inositol_phosphate_patterns:
        if mol.HasSubstructMatch(pattern):
            has_inositol_phosphate = True
            break

    if not has_inositol_phosphate:
        return False, "Inositol ring lacks phosphate group or incorrect substitution"



    # Check for ester groups (-O-C(=O)-) and fatty acid chains
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
      return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains. Using a looser criteria
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
       return False, "Chains too short to be fatty acids"

    # Molecular weight - typical range is > 500 - using a lower threshold
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
       return False, "Molecular weight too low for PIP"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 20:
       return False, "Too few carbons"
    if p_count < 1:
      return False, "Must have at least one phosphorus"
    if o_count < 10:
      return False, "Too few oxygen atoms"


    return True, "Contains a glycerol backbone, a phosphate, an inositol ring with at least one phosphate, and fatty acid chains"