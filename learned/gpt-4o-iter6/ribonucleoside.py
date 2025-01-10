"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is defined by the presence of a nucleobase attached to a D-ribose sugar.

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

    # Ribose sugar pattern adjusted for more flexibility (account for D-ribose with common stereocenters)
    ribose_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H](O1)CO)O)O)O")
    
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar ring found"

    # Nucleobase patterns recognizing common purines and pyrimidines
    # Here, care is taken to accommodate various substitutions and tautomers
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1nc[nH]c2"),  # Purine
        Chem.MolFromSmarts("c1ncnc2[nH]cnc12"),  # Pyrimidine and variations
        Chem.MolFromSmarts("c1nc[nH]c2[nH]cnc12"),  # Extended purine match
        Chem.MolFromSmarts("c1c[nH]cnc1"),  # Additional heterocyclic components
    ]

    # Check if the molecule has any nucleobase attached to ribose
    nucleobase_found = any(mol.HasSubstructMatch(nb_pattern) for nb_pattern in nucleobase_patterns)
    if not nucleobase_found:
        return False, "No nucleobase found attached to ribose"

    return True, "Contains ribose sugar with an attached nucleobase"

# Example testing of the function
smiles_examples = [
    "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cc(CNCC(O)=O)c(=O)[nH]c1=O",  # 5-(carboxymethylaminomethyl)uridine
    "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",  # guanosine
    "C[C@H](CNC(=O)CC[C@]1(C)[C@@H](CC(N)=O)[C@H]2N3C1=C(C)C1=[N+]4C(=CC5=[N+]6C(=C(C)C7=[N+]([C@]2(C)[C@@](C)(CC(N)=O)[C@@H]7CCC(N)=O)[Co--]346C[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@](C)(CC(N)=O)[C@@H]5CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O)OP(O)(O)=O"  # adenosylcobinamide phosphate
]

for smiles in smiles_examples:
    is_ribo, reason = is_ribonucleoside(smiles)
    print(f"SMILES: {smiles}, is_ribonucleoside: {is_ribo}, Reason: {reason}")