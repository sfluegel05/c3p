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
    thiamine_pattern = Chem.MolFromSmarts("C1=CSC2=[N+]1CCO")  # captures thiazole and pyrimidine parts
    riboflavin_pattern = Chem.MolFromSmarts("C1=NC2=C(C=C1C)N=C(NC2=O)NC3=C(N=C(N=CC3=O)N(C)C)C")  # extends pattern for isoalloxazine
    niacin_pattern = Chem.MolFromSmarts("n1ccc(C(=O)[O,NH1])cn1")  # niacin or nicotinamide subtype
    pantothenic_acid_pattern = Chem.MolFromSmarts("OC(=O)C(C(C(=O)O)NC(=O)C(C)O)CO")  # includes pantolactone
    pyridoxine_pattern = Chem.MolFromSmarts("C1=NC=C(C(=C1CO)CO)CO")  # pyridine with hydroxymethyl groups
    biotin_pattern = Chem.MolFromSmarts("C1(C(=O)NC2C(S1)CCC2)CCC(=O)[O,N]C")  # sulfur and ureido structures
    folate_pattern = Chem.MolFromSmarts("NC1=NC2=[NH]c3cnc(N)c3n2nc1=O")  # captures pteridine core
    cobalamin_pattern = Chem.MolFromSmarts("[Co]1N2C=CC=C2C=C1")  # simplified to catch corrin ring motifs and cobalt

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