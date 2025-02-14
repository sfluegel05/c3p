"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_6_aminopurine(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine based on its SMILES string.
    A 6-aminopurine is any compound having 6-aminopurine (adenine) as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define 6-aminopurine (adenine) pattern
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    
    # Check if molecule contains adenine pattern
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) moiety"
    else:
        return False, "Does not contain 6-aminopurine (adenine) moiety"

# Example usage
smiles1 = "Nc1ncnc2n(cnc12)C3OC(COP(=O)(O)OP(=O)(O)OCC4OC(n5cnc6c(N)ncnc56)C(O)C4O)C(O)C3O"  # Acetyl-CoA
smiles2 = "CCc1ccccc1"  # Ethylbenzene
smiles3 = "Nc1ncnc2[nH]cnc12"  # 6-methyladenine
smiles4 = "OC(CN1C2=NC=NC(N)=C2N=C1)C(O)C(O)=O"  # D-erythro-Eritadenine

print(is_6_aminopurine(smiles1))  # Output: (True, 'Contains 6-aminopurine (adenine) moiety')
print(is_6_aminopurine(smiles2))  # Output: (False, 'Does not contain 6-aminopurine (adenine) moiety')
print(is_6_aminopurine(smiles3))  # Output: (True, 'Contains 6-aminopurine (adenine) moiety')
print(is_6_aminopurine(smiles4))  # Output: (True, 'Contains 6-aminopurine (adenine) moiety')