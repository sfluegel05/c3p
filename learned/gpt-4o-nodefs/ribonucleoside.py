"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Pattern for ribose sugar - 5-membered ring with C-C-C-C-O and three OH groups
    ribose_pattern = Chem.MolFromSmarts("C1OC([H,C])([H,C])C(O)C1O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar moiety found"
    
    # Pattern for glycosidic bond to base (simplified) - base nitrogen bonded to ribose
    base_linkage_pattern = Chem.MolFromSmarts("[nH]1[c,n,o][c,n,o][nH][c,n,o]1[C@H]1O[C@H]([C@@H](O)[C@H]1O)C")
    if not mol.HasSubstructMatch(base_linkage_pattern):
        return False, "No base linked to ribose sugar"
    
    return True, "Contains ribose sugar with glycosidic linkage to a base"

# Example usage
smiles_examples = [
    "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",  # guanosine
    "CCCCCOC(=O)Nc1nc(=O)n(cc1F)[C@@H]1O[C@H](C)[C@@H](O)[C@H]1O",  # capecitabine
    "CSC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12"  # S-adenosyl-3-thiopropylamine
]

for smiles in smiles_examples:
    is_ribo, reason = is_ribonucleoside(smiles)
    print(f"SMILES: {smiles}\nIs Ribonucleoside: {is_ribo}, Reason: {reason}\n")