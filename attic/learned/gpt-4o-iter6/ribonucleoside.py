"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is identified by the presence of a nucleobase attached to a D-ribose sugar.

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

    # Define ribose ring pattern (D-ribose)
    ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@@H](CO)[C@H]1O)")

    # Check for ribose sugar
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar ring found"

    # Define common nucleobase backbone patterns
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1nc[nH]c2"),  # purine (adenine, guanine)
        Chem.MolFromSmarts("c1ncnc2[nH]cnc12"),  # pyrimidine (cytosine, uracil)
    ]

    # Check for nucleobase attachment
    nucleobase_found = any(mol.HasSubstructMatch(nb_pattern) for nb_pattern in nucleobase_patterns)
    if not nucleobase_found:
        return False, "No nucleobase found attached to ribose"

    return True, "Contains ribose sugar with attached nucleobase"

# Example testing of the function
smiles_examples = [
    "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cc(CNCC(O)=O)c(=O)[nH]c1=O",  # 5-(carboxymethylaminomethyl)uridine
    "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"  # guanosine
]

for smiles in smiles_examples:
    is_ribo, reason = is_ribonucleoside(smiles)
    print(f"SMILES: {smiles}, is_ribonucleoside: {is_ribo}, Reason: {reason}")