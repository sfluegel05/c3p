"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose or deoxyribose sugar linked to a purine or pyrimidine base,
    with phosphorylation at the C-5 position of the sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a ribosyl or deoxyribosyl sugar with flexibility in stereochemistry
    # and compatibility with nucleoside linking.
    sugar_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@H](O)[C@H](CO)O1")  # Accept natural variations
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Define a broader pattern for phosphate group/s attached to the 5' carbon of the sugar.
    phosphate_pattern = Chem.MolFromSmarts("O[C@H]1COP(O)(=O)O[*1]")  # Allow for extended phosphate groups
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found at the 5' position of the sugar"

    # Check for presence of purine or pyrimidine base, assume it's linked to the sugar
    # Accepting purine/pyrimidine base connected through N-glycosidic bond
    base_pattern = Chem.MolFromSmarts("n1cnc2c1[nH]c(=O)[nH]c2=O")  # Simplified pyrimidine (e.g. uracil)
    purine_pattern = Chem.MolFromSmarts("N1C=NC2=C1N=CN=C2")  # Simplified purine (e.g. adenine)
    if not mol.HasSubstructMatch(base_pattern) and not mol.HasSubstructMatch(purine_pattern):
        return False, "No nucleobase found linked to sugar"

    return True, "Contains ribosyl or deoxyribosyl sugar, phosphate group at the 5' position, and nucleobase"