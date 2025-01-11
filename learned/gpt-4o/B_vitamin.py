"""
Classifies: CHEBI:75769 B vitamin
"""
from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins include Thiamine, Riboflavin, Niacin, Pantothenic acid,
    Pyridoxine, Biotin, Folic acid, and Cobalamin.

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

    # SMARTS patterns for key B vitamins:
    patterns = {
        "B1": Chem.MolFromSmarts("C[C@H]1CN(C2=CN=C(C)N=C2N)C=C1"), # Thiamine
        "B2": Chem.MolFromSmarts("C1(=CC=C2C(O)=CN=C(NC)C2=CN=1)"), # Riboflavin
        "B3": Chem.MolFromSmarts("OC(=O)c1cccnc1"), # Nicotinic acid
        "B5": Chem.MolFromSmarts("CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O"), # Pantothenic acid
        "B6": Chem.MolFromSmarts("CC1=NC=C(C=C1O)O"), # Pyridoxine
        "B7": Chem.MolFromSmarts("NC1N(A)S(A)OC1"), # Biotin (simplified key structure)
        "B9": Chem.MolFromSmarts("Nc1nc2NCCNC2C(=O)N1"), # Folic acid
        "B12": Chem.MolFromSmarts("CNC1=CN=C(C=C1)C"), # Cobalamin, indicative corrin structure
    }

    # Check for substructures
    for vitamin, pattern in patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f"Matches pattern of vitamin {vitamin}"
    
    return False, "Does not match any known B vitamin patterns"