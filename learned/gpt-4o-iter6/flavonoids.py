"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extended SMARTS patterns for flavonoids covering diverse structures
    flavonoid_smarts = [
        "c1c(O)cc(O)c2c1oc3c(c2=O)ccc(=O)c3",  # General flavonoid backbone
        "c1cc(O)c2c(c1)oc3c(c2=O)cc(O)c(c3)O",  # Flavonol backbone
        "c1cc2c(oc=c(O)c2=O)c3ccc(O)cc13",  # Isoflavonoid structure
        "c1(C=O)cc(O)cc1-c2ccoc2",  # Chalcone structure
        "c1c(O)cc(O)c2c1oc(c(c2=O)C)C",  # Flavanone structure
        "[OH]c1c(O)cc(C(=O)C=C2C=CC=CC=2)ccc1",  # Phenylpropane structure
        "c1cc2c(c(c1O)O)oc(=O)c3c2c(ccc3)O",  # Flavone structure
        "c1c(O)cc2c(c1)C(=O)c3c(O)cc(O)cc3O2",  # Biflavonoid structure
        "c1cc(O)c2c(c1)C(=O)c3c(O)cc(O)cc3O2",  # Naringenin-like
        "[O][C@H]1COC(c2cc(O)cc(O)c2)c3cc(O)cc(O)c13"  # Glycosylated flavonoid
    ]

    # Check for matches with any of the patterns
    for pattern in flavonoid_smarts:
        flavonoid_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(flavonoid_pattern):
            return True, f"Contains a recognized flavonoid structure: {pattern}"

    return False, "Does not contain a recognized flavonoid structure"

# Example usage:
smiles_examples = [
    "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC(O)=C(O)C(O)=C4)CO",
    "C=1C(=CC2=C(C1O)C([C@@H]([C@H](O2)C=3C=CC(=CC3)O)O[C@H]4[C@@H]([C@@H]([C@H]([C@@H](O4)CO)O)O)O)=O)O",
    "C1(C2=C(O[C@@H](C1)C3=C(C=C(C=C3)O)O)C(=C(C=C2O)O)C[C@@H](CCC(=C)C)C(=C)C)=O"
]

for smi in smiles_examples:
    result, reason = is_flavonoids(smi)
    print(f"SMILES: {smi}\nIs Flavonoid: {result}, Reason: {reason}\n")