"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    Since 'semisynthetic derivative' involves context about origin, this function uses some heuristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule may be a semisynthetic derivative, False otherwise
        str: Reason for classification or note on heuristic limitations
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic: Check for presence of unique functional groups often introduced synthetically
    # Such as esters, ethers, or phosphate groups which may indicate synthetic modification
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")

    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    # Count number of functional groups
    num_functional_groups = 0
    if has_ester:
        num_functional_groups += 1
    if has_ether:
        num_functional_groups += 1
    if has_phosphate:
        num_functional_groups += 1

    # If two or more types of synthetic lips are present, possibly semisynthetic
    if num_functional_groups >= 2:
        return True, "Multiple synthetic-looking groups detected"

    # Without explicit synthesis pathway info, classification is highly speculative
    return False, "Unable to determine synthesis origin from SMILES"

# Examples to test the function
print(is_semisynthetic_derivative("CO[C@H]1C[C@H](O[C@H]2[C@H](C)[C@@H](O[C@@H]3O[C@H](C)C[C@@H]([C@H]3OC(C)=O)N(C)C)[C@@H](C)C[C@@]3(CO3)C(=O)[C@H](C)[C@@H](OC(C)=O)[C@@H](C)[C@@H](C)OC(=O)[C@@H]2C)O[C@@H](C)[C@@H]1OC(C)=O"))  # Troleandomycin
print(is_semisynthetic_derivative("CCCCn1c2cc(OC)ccc2c2ccnc(C)c12"))  # N-butylharmine