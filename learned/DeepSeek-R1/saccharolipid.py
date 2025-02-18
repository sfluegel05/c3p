"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI: saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    Saccharolipids contain carbohydrate (sugar) and lipid components covalently linked via ester/amide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced carbohydrate detection: Match common pyranose/furanose structures with multiple hydroxyls
    # Modified to include common sugar substitutions (e.g. amino sugars, deoxy)
    sugar_core = Chem.MolFromSmarts("""
        [C;r5,r6]1@[C;r5,r6]@[C;r5,r6]@[C;r5,r6]@[C;r5,r6][OH1,OH0,C(=O),O;r5,r6]1
        """)
    if not mol.HasSubstructMatch(sugar_core):
        return False, "No carbohydrate core detected (pyranose/furanose structure)"

    # Lipid detection: Long chain (>=8 carbons) connected via ester/amide to sugar
    # SMARTS pattern for ester/amide connected to both sugar oxygen and long chain
    lipid_linker = Chem.MolFromSmarts("""
        [O;R][C](=O)[C][C][C][C][C][C][C,C] 
        |^1: ester/amide connected to sugar oxygen,
        chain_start: [C]=O connected to at least 7 more carbons
        """)
    if not mol.HasSubstructMatch(lipid_linker):
        # Check for amide linkages common in lipid A structures
        lipid_amide = Chem.MolFromSmarts("""
            [C](=O)[N][C][C][C][C][C][C][C,C][C,C] 
            |^1: amide connected to long chain
            """)
        if not mol.HasSubstructMatch(lipid_amide):
            return False, "No lipid chain (>=8 carbons) connected via ester/amide"

    # Exclude phosphate-linked structures (common in phospholipids)
    phosphate = Chem.MolFromSmarts("[O]P(=O)([O])[O]")
    if mol.HasSubstructMatch(phosphate):
        # Allow phosphate only if not between sugar and lipid
        # Check if phosphate is part of the linker between sugar and lipid
        ester_amide_with_phos = Chem.MolFromSmarts("[O;R]-P(=O)([O])-[O]-C(=O)")
        if mol.HasSubstructMatch(ester_amide_with_phos):
            return False, "Phosphate bridge between sugar and lipid"

    return True, "Contains carbohydrate linked to lipid via ester/amide"