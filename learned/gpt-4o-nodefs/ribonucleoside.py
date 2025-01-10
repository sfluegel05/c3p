"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside has a base (purine or pyrimidine) attached to a ribose sugar.

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

    # General pattern for ribose - 5-membered carbon ring with multiple hydroxyl or oxygen attachments
    ribose_pattern = Chem.MolFromSmarts("C1OC[C@H](O)[C@@H]1O")
    deoxyribose_pattern = Chem.MolFromSmarts("C1OC[C@H](O)C1")

    ribose_match = mol.HasSubstructMatch(ribose_pattern)
    deoxyribose_match = mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not ribose_match and not deoxyribose_match:
        return False, "No ribose or deoxyribose sugar moiety found"

    # Allow variance in bases
    purine_pattern = Chem.MolFromSmarts("n1cnc2[nH]cnc12")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncn[C@H]1")

    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No purine or pyrimidine moiety found"

    # Ensure there is a glycosidic bond to the sugar
    linkage_pattern = Chem.MolFromSmarts("[n][C@H]1O")
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No glycosidic bond connecting base and sugar"

    return True, "Contains ribose or deoxyribose sugar with glycosidic linkage to a purine or pyrimidine base"

# Example usage
smiles_examples = [
    "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",  # guanosine
    "CCCCCOC(=O)Nc1nc(=O)n(cc1F)[C@@H]1O[C@H](C)[C@@H](O)[C@H]1O",  # capecitabine
    "CSC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"  # S-adenosyl-3-thiopropylamine
]

for smiles in smiles_examples:
    is_ribo, reason = is_ribonucleoside(smiles)
    print(f"SMILES: {smiles}\nIs Ribonucleoside: {is_ribo}, Reason: {reason}\n")