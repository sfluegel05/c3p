"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains the -C#N group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyano group pattern (C#N)
    cyano_pattern = Chem.MolFromSmarts("C#N")
    if mol.HasSubstructMatch(cyano_pattern):
        return True, "Contains a cyano group (C#N)"
    else:
        return False, "No cyano group found"

# Example usage
smiles_list = [
    "COC1=CC=CC=C1NC(=O)N2CCCCN3[C@H](C2)[C@@H]([C@H]3CO)C4=CC=C(C=C4)C5=CC=CC(=C5)C#N",
    "CCC(C)C#N",
    "ClC(Br)C#N",
    "CN(C)CCC[C@@]1(OCc2cc(ccc12)C#N)c1ccc(F)cc1"
]

for smiles in smiles_list:
    result, reason = is_nitrile(smiles)
    print(f"SMILES: {smiles} -> Is Nitrile: {result}, Reason: {reason}")