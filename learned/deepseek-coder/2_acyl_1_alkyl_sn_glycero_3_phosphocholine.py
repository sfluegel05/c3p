"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:73005 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    The molecule must have a glycerol backbone with a phosphocholine group at the 3-position,
    an alkyl group at the 1-position, and an acyl group at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[NX4+](C)(C)CCOP([O-])(=O)")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for the glycerol backbone with the correct stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COP([O-])(=O)OCC[N+](C)(C)C)(COC)OC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-glycero backbone with phosphocholine at 3-position found"

    # Check for the presence of an alkyl group at the 1-position
    # The pattern now allows for any length of alkyl chain
    alkyl_pattern = Chem.MolFromSmarts("[CX4][OX2][C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No alkyl group at 1-position found"

    # Check for the presence of an acyl group at the 2-position
    # The pattern now allows for any length of acyl chain
    acyl_pattern = Chem.MolFromSmarts("[CX4][OX2][C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)[CX4]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group at 2-position found"

    # Check molecular weight - these compounds typically have a high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"
    if o_count < 6:
        return False, "Too few oxygens for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"

    return True, "Contains sn-glycero backbone with phosphocholine at 3-position, alkyl at 1-position, and acyl at 2-position"