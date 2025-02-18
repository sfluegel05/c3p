"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: CHEBI:17515 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol consists of a glycerol backbone with two fatty acid esters
    and a phosphatidyl group (glycerol connected via phosphate).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of phosphorus (phosphate group)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphate group found"

    # Check for at least two ester groups (O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Core structure: Phosphate connecting two glycerol units
    # Pattern matches phosphate (P) connected to two oxygen atoms that are part of glycerol backbones
    core_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-OP(=O)([O-])-[O][CH2]-[CH](-[CH2])-[CH2])-[CH2]"
    )
    if not mol.HasSubstructMatch(core_pattern):
        # Try alternative pattern accounting for possible stereochemistry and connectivity
        core_pattern = Chem.MolFromSmarts(
            "[CH2](-[O]-P(=O)([O-])-[O]-[CH2]-[CH](-[CH2])-[CH2])-[CH](-[CH2])-[CH2]"
        )
        if not mol.HasSubstructMatch(core_pattern):
            return False, "Glycerol-phosphate-glycerol backbone not found"

    # Verify at least two hydroxyl groups in the head glycerol (non-ester oxygens)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 2:
        return False, "Insufficient hydroxyl groups in structure"

    # Check for fatty acid chains (minimum 8 carbons each)
    # Approximated by checking rotatable bonds > 10 (indicative of long chains)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 10:
        return False, f"Only {rot_bonds} rotatable bonds, insufficient for fatty acid chains"

    # Molecular weight check (typical PG > 600 but some examples may be lower)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Contains glycerol-phosphate-glycerol backbone with two fatty acid esters"