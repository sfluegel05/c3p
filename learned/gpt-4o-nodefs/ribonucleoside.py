"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside has a nucleobase (purine or pyrimidine) attached to a ribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for ribose-like sugar
    ribose_pattern = Chem.MolFromSmarts("C1OC(CO)[C@@H](O)[C@H]1O")
    
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar moiety found"

    # Patterns for nucleobases
    purine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ccncn1")

    # Check for nucleobase
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No purine or pyrimidine nucleobase found"

    # Verify glycosidic bond, i.e., checks linkage between sugar and base
    linkage_pattern = Chem.MolFromSmarts("c1[nH]c2c(ncnc12)[C@H]3O")
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No proper glycosidic linkage between base and sugar"

    return True, "Contains ribose sugar with glycosidic linkage to a purine or pyrimidine base"

# Example usage
smiles_examples = [
    "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",  # guanosine
    "CCCCCOC(=O)Nc1nc(=O)n(cc1F)[C@@H]1O[C@H](C)[C@@H](O)[C@H]1O",  # capecitabine
    "CSC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"  # S-adenosyl-3-thiopropylamine
]

for smiles in smiles_examples:
    is_ribo, reason = is_ribonucleoside(smiles)
    print(f"SMILES: {smiles}\nIs Ribonucleoside: {is_ribo}, Reason: {reason}\n")