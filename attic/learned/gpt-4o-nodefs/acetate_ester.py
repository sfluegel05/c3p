"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is characterized by the presence of an acetate group (CH3COO-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetate ester pattern (C(=O)OC with a methyl group)
    acetate_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    methyl_group_pattern = Chem.MolFromSmarts("C([CH3])") # Specific methyl group check

    # Find the ester linkage
    if not mol.HasSubstructMatch(acetate_ester_pattern):
        return False, "No acetate ester linkage found"

    # Confirm the ester linkage includes the specific methyl group
    if not mol.HasSubstructMatch(methyl_group_pattern):
        return False, "Ester linkage does not include methyl group"

    return True, "Molecule contains an acetate ester group"

# Examples for testing
example_smiles = [
    "CC(=O)OCC",  # Simple acetate ester
    "CC(=O)Oc1ccccc1",  # Phenyl acetate
    "C1=CC=C(C=C1)OCC(=O)C",  # Does not match because of incorrect arrangement
    "C1=CC(=CC=C1C(=O)O)N",  # Incorrect - not an acetate ester
]

for smiles in example_smiles:
    is_acetate, reason = is_acetate_ester(smiles)
    print(f"SMILES: {smiles}\nIs Acetate Ester: {is_acetate}\nReason: {reason}\n")