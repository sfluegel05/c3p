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
        bool: True if the molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key structural features of B vitamins
    try:
        thiamine_pattern = Chem.MolFromSmarts("Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1")  # Captures thiazole and pyrimidine rings
        riboflavin_pattern = Chem.MolFromSmarts("Cc1cc2nc3c(nc(=O)[n-]c3=O)n(CC[C@@H](O)[C@H](O)[C@H](O)CO)c2cc1C")  # Isoalloxazine core
        niacin_pattern = Chem.MolFromSmarts("n1ccc(C(=O)[O,N])cn1")  # Nicotinic acid (niacin)
        pantothenic_acid_pattern = Chem.MolFromSmarts("OC(=O)C(C(C(=O)O)NC(=O)C(C)O)CO")  # Pantothenic acid
        pyridoxine_pattern = Chem.MolFromSmarts("C1=NC=C(C(=C1CO)CO)CO")  # Pyridoxine structure
        biotin_pattern = Chem.MolFromSmarts("C1(C(=O)NC2C(S1)CCC2)CCC(=O)[O,N]")  # Biotin structural elements
        folate_pattern = Chem.MolFromSmarts("Nc1nc2NCC(CNc3ccc(cc3)C(=O)NC[C@@H](CCC(O)=O)C(O)=O)Nc2c(=O)[nH]1")  # Simplified folate structure
        # Note: Cobalamin (B12) is complex and involves metal coordination bonds which may not be easily captured by simple SMARTS

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

        # Cobalamin check could be added with a specialized structure pattern from a library if needed.

        return False, "No match for B vitamin structures found in SMILES"
    except Exception as e:
        return None, f"Error during pattern matching: {str(e)}"