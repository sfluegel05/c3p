"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    A 1-phosphatidyl-1D-myo-inositol is a phosphatidylinositol where the inositol moiety is the 1D-myo isomer,
    and the phosphatidyl group is located at position 1. It should have two fatty acid chains attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 1D-myo-inositol backbone
    inositol_pattern = Chem.MolFromSmarts("[C@H]1(O[C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)O")
    inositol_match = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_match:
        return False, "No 1D-myo-inositol backbone found"

    # Check for phosphate group at position 1
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)O")
    phosphate_match = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_match:
        return False, "No phosphate group at position 1 found"

    # Check for two fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern, maxMatches=2)
    if len(fatty_acid_matches) != 2:
        return False, f"Found {len(fatty_acid_matches)} fatty acid chains, expected exactly 2"

    # Check for long fatty acid chains (>= 12 carbon atoms)
    fatty_acid_chains = []
    for match in fatty_acid_matches:
        chain = []
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain.append(atom)
        if len(chain) >= 12:
            fatty_acid_chains.append(chain)
    if len(fatty_acid_chains) != 2:
        return False, "Fatty acid chains are too short"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 25 or o_count < 9:
        return False, "Insufficient carbon and oxygen atoms for a 1-phosphatidyl-1D-myo-inositol"

    return True, "Molecule is a 1-phosphatidyl-1D-myo-inositol with two fatty acid chains attached"