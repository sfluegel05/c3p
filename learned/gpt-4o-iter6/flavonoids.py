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
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flavonoids are often characterized by their phenyl-benzopyranone structures:
    flavonoid_smarts = [
        "c1cc2oc(c(c2c(c1)O)O)-c3ccccc3",  # Flavonoids base pattern (flavanone)
        "c1cc2c(cc1)c(=O)c3c(o2)cc(O)c(O)c3",  # Isoflavonoid pattern (isoquinolinone)
        "c1cc(O)c2c(c1)oc3c(c2=O)cc(O)c(c3)O"  # Neoflavonoid pattern
    ]
    
    for pattern in flavonoid_smarts:
        flavonoid_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(flavonoid_pattern):
            return True, f"Contains a recognized flavonoid backbone: {pattern}"

    return False, "Does not contain a recognized flavonoid backbone"

# Example usage:
smiles_examples = [
    "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC=2C=C3C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO)=CC(O)=CC3=[O+]C2C5=CC=C(O)C=C5)CO",
    # Add more SMILES strings for testing
]

for smi in smiles_examples:
    result, reason = is_flavonoids(smi)
    print(f"SMILES: {smi}\nIs Flavonoid: {result}, Reason: {reason}\n")