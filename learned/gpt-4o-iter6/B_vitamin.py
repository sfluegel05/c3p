"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamins
"""

from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key structural features of B vitamins
    thiamine_pattern = Chem.MolFromSmarts("C1=CSC([N+]=C1)CCCOP")
    riboflavin_pattern = Chem.MolFromSmarts("C1=C2C=CC([N]=C3[C@H]([C@H](O)[C@@H](O)[C@H](O)CO)NC3=O)=N2N=C1")
    niacin_pattern = Chem.MolFromSmarts("c1cccnc1C(=O)O")
    pantothenic_acid_pattern = Chem.MolFromSmarts("C(C(CCC(=O)O)NC(=O)C(C)C)(O)O")
    pyridoxine_pattern = Chem.MolFromSmarts("C1=NC=C(C(=C1CO)CO)CO")
    biotin_pattern = Chem.MolFromSmarts("C1(C(=O)NC2C(S1)CCC2)CCC(C(=O)O)O")
    folate_pattern = Chem.MolFromSmarts("c1nc(NCC2=CC=C(C=C2)C(=O)NC(CC(=O)O)C(=O)O)[nH]c3[nH]cnc1c3=O")
    cobalamin_pattern = Chem.MolFromSmarts("[H][C@](C)(P(=O)(O)OCC[C@@H](CO)[C@H](O[C@@H]1O[C@H](n2cnc3c(ncnc23)N)[C@H]([C@H]1O)O)O)c1ccc(cc1)C#N")

    # Check against each B vitamin pattern
    if mol.HasSubstructMatch(thiamine_pattern):
        return True, "Matches thiamine (B1) structure"
    if mol.HasSubstructMatch(riboflavin_pattern):
        return True, "Matches riboflavin (B2) structure"
    if mol.HasSubstructMatch(niacin_pattern):
        return True, "Matches niacin (B3) structure"
    if mol.HasSubstructMatch(pantothenic_acid_pattern):
        return True, "Matches pantothenic acid (B5) structure"
    if mol.HasSubstructMatch(pyridoxine_pattern):
        return True, "Matches pyridoxine (B6) structure"
    if mol.HasSubstructMatch(biotin_pattern):
        return True, "Matches biotin (B7) structure"
    if mol.HasSubstructMatch(folate_pattern):
        return True, "Matches folate (B9) structure"
    if mol.HasSubstructMatch(cobalamin_pattern):
        return True, "Matches cobalamin (B12) structure"

    return False, "No match for B vitamin structures found in SMILES"