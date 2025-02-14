"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate has a glycerol backbone, a phosphate group,
    a myo-inositol ring, and fatty acid chains. Crucially it contains an additional phosphate group on the inositol ring.

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

    # 1. Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # 2. Check for phosphate group connected to the glycerol backbone via an oxygen.
    #    This phosphate must also be an ester
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[PX4](=[OX1])(-[OX2])")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No phosphate group attached to glycerol backbone"


    # 3. Check for myo-inositol ring (six-membered ring with 6 hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("C1[C@@H]([C@@H]([C@@H]([C@H]([C@@H]1O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
       return False, "No inositol ring found"

    # 4. Check for phosphate group(s) on the inositol ring. Use a flexible pattern.
    inositol_phosphate_pattern1 = Chem.MolFromSmarts("C1[C@H]([C@H]([C@H]([C@@H]([C@H]1O)O)O)[OX2]-[PX4](=[OX1])([OX2,O-])[OX2,O-])") #direct
    inositol_phosphate_pattern2 = Chem.MolFromSmarts("C1[C@H]([C@H]([C@H]([C@@H]([C@H]1O)[OX2]-[#6])O)O)[OX2]-[PX4](=[OX1])([OX2,O-])[OX2,O-]") # single atom between ring and P
    if not (mol.HasSubstructMatch(inositol_phosphate_pattern1) or mol.HasSubstructMatch(inositol_phosphate_pattern2)):
        return False, "Inositol ring lacks phosphate group"

    # 5. Check that the glycerol phosphate is linked to the inositol ring
    glycerol_inositol_linkage_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[PX4](=[OX1])(-[OX2])-[OX2]-[#6]1[#6]([#6]([#6]([#6]([#6]1[OX2])OX2)OX2)OX2)OX2")
    if not mol.HasSubstructMatch(glycerol_inositol_linkage_pattern):
        return False, "Glycerol phosphate not linked to inositol"

    # 6. Check for ester groups (-O-C(=O)-) and fatty acid chains (at least 2)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2 :
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"


    # Check for a minimum number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
      return False, "Chains too short to be fatty acids"
      
    # Check for a minimum molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
      return False, "Molecular weight too low for PIP"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 20:
       return False, "Too few carbons for a PIP"
    if p_count < 2:
      return False, "Must have at least two phosphorus for a PIP"
    if o_count < 10:
      return False, "Too few oxygen atoms"

    return True, "Contains a glycerol backbone, a phosphate, an inositol ring with at least one phosphate, and fatty acid chains"