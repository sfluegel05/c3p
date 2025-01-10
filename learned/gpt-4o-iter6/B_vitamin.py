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
    thiamine_pattern = Chem.MolFromSmarts("C1=CSC([N+]=C1)CCO")  # Adjusted to capture phosphorylated derivatives
    riboflavin_pattern = Chem.MolFromSmarts("C1=NC2=C(C=C1C)N=C(NC2=O)N")  # Adjusted for simpler ring matches
    niacin_pattern = Chem.MolFromSmarts("c1cccnc1C(=O)[O,N]")  # Can match broader functional groups
    pantothenic_acid_pattern = Chem.MolFromSmarts("C(C(C(C(=O)O)NC(=O)C(C)C)O)O")  # Retained initially
    pyridoxine_pattern = Chem.MolFromSmarts("C1=NC=C(C(=C1CO)CO)CO")  # Retained initially
    biotin_pattern = Chem.MolFromSmarts("C1(C(=O)NC2C(S1)CCC2)CCC(C(=O)O)O")  # Revalidate for ring closure
    folate_pattern = Chem.MolFromSmarts("c1nc(NC2=CC=CC=C2)nc(N)c1=O")  # Simplified p-aminobenzoate portion
    cobalamin_pattern = Chem.MolFromSmarts("C[cH]1cC[nH]c2c[nH][nH]c21")  # Simplify and expand for B12 complexity

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