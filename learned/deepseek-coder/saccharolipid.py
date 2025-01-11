"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
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

    # Check for lipid-like features (long hydrocarbon chains, ester bonds)
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found (lipid-like feature missing)"

    # Check for carbohydrate-like features using more flexible patterns
    # Look for any 5- or 6-membered rings with multiple oxygens
    sugar_pattern = Chem.MolFromSmarts("[C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][O]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 1:
        sugar_pattern = Chem.MolFromSmarts("[C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][O]")
        sugar_matches = mol.GetSubstructMatches(sugar_pattern)
        if len(sugar_matches) < 1:
            return False, "No sugar rings found (carbohydrate-like feature missing)"

    # Check for glycosidic bonds (C-O-C between sugar and any other part)
    glycosidic_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C;H1,H2]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic bonds found (carbohydrate-like feature missing)"

    # Check for long aliphatic chains (lipid-like feature)
    aliphatic_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if len(aliphatic_matches) < 2:
        return False, "Not enough long aliphatic chains found (lipid-like feature missing)"

    # Additional check for connectivity between lipid and carbohydrate parts
    # Look for ester bonds connected to carbohydrate rings
    lipid_carb_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][C;H1,H2][C;H1,H2][C;H1,H2][C;H1,H2][O]")
    lipid_carb_matches = mol.GetSubstructMatches(lipid_carb_pattern)
    if len(lipid_carb_matches) < 1:
        return False, "No connection between lipid and carbohydrate parts"

    # Check molecular weight - saccharolipids typically >500 Da
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for saccharolipid"

    # Check number of rotatable bonds to ensure sufficient flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds for saccharolipid"

    return True, "Contains both lipid-like and carbohydrate-like features with proper connectivity"