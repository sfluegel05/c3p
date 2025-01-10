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

    # Improved SMARTS patterns for key B vitamins
    patterns = {
        "B1": Chem.MolFromSmarts("C[C@H]1CN(C2=CN=C(C)N=C2N)C=C1"),  # Thiamine, emphasize structure with the thiazole ring and aminopyrimidine
        "B2": Chem.MolFromSmarts("NC1=NC2=C(NC=N2)C3=C(C(O)=CN=C13)"),  # Riboflavin, capture the isoalloxazine ring
        "B3": Chem.MolFromSmarts("C1=NC=CC=C1C(=O)O"),  # Nicotinic acid, capture the pyridine ring
        "B5": Chem.MolFromSmarts("CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)O"),  # Pantothenic acid
        "B6": Chem.MolFromSmarts("CC1=NC=C(CN)C(O)=C1"),  # Pyridoxine, improved capture of the pyridine ring with hydroxyl
        "B7": Chem.MolFromSmarts("NC1CS[C@H]2N1C=O"),  # Biotin, more detailed to include urea and tetrahydrothiophene
        "B9": Chem.MolFromSmarts("Nc1nc2ccc(O)c(=O)n2[nH]1"),  # Folic acid, focus on the pteridine and p-aminobenzoic acid portions
        "B12": Chem.MolFromSmarts("CNC([Co])C1=NC=CC(=N1)C"),  # Cobalamin, focus on large corrin with cobalt 
    }

    # Check for substructures
    for vitamin, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches pattern of vitamin {vitamin}"
    
    return False, "Does not match any known B vitamin patterns"