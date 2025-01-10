"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Look for flavonoid backbone pattern: C6-C3-C6
    # Using a simplified structure of [Ar]c1cc([Ar])c([Ar])cc1 where [Ar] indicates aromaticity
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(O)c2c(c1)oc(-c1ccccc1)c(O)c2=O")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Does not contain the core C6-C3-C6 flavonoid backbone"

    # Check for other common flavonoid features (optional additional validation)
    # We won't go into full subclasses detail here due to complexity.

    return True, "Contains core C6-C3-C6 flavonoid backbone"

# Example usage:
smiles_examples = [
    "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC=2C=C3C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO)=CC(O)=CC3=[O+]C2C5=CC=C(O)C=C5)CO"
    # add more SMILES strings for testing
]

for smi in smiles_examples:
    result, reason = is_flavonoids(smi)
    print(f"SMILES: {smi}\nIs Flavonoid: {result}, Reason: {reason}\n")