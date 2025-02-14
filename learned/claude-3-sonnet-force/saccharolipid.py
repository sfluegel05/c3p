"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI:26094 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid that contains a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbohydrate patterns
    sugar_pattern = Chem.MolFromSmarts("[OX2H][CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[OX2H,CX4H]")  # cyclic structure with multiple hydroxyl groups
    sugar_matches = mol.HasSubstructMatch(sugar_pattern)
    if not sugar_matches:
        return False, "No carbohydrate moiety found"

    # Look for lipid patterns
    lipid_patterns = [Chem.MolFromSmarts("[CX3](=O)[OX2H]"),  # ester
                      Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H2][CX4H2][CX4H2]"),  # alkyl chain
                      Chem.MolFromSmarts("[SX4+2]([OX1-])([OX1-])[OX2H]"),  # sulfate ester
                      Chem.MolFromSmarts("[CX4H2]=[CX3H]([CX4H2])[CX4H2]"),  # prenyl group
                      Chem.MolFromSmarts("[NX3H2]")]  # amine group
    lipid_matches = any(mol.HasSubstructMatch(pattern) for pattern in lipid_patterns)
    if not lipid_matches:
        return False, "No lipid moiety found"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for saccharolipid"

    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if o_count < 5 or c_count < 10:
        return False, "Insufficient oxygen or carbon atoms for saccharolipid"

    return True, "Contains both carbohydrate and lipid moieties"