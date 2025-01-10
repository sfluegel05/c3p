"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a nucleotide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify at least one phosphate group, considering variations (mono, di, tri, cyclic)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 1:
        return False, "No phosphate group found"
    
    # Check for the presence of a ribose or deoxyribose sugar ring
    ribose_pattern = Chem.MolFromSmarts("C[C@H](O)[C@H](O)[C@H](CO)O")
    deoxyribose_pattern = Chem.MolFromSmarts("C[C@H](O)[C@H](O)[C@H](CO)")

    if not mol.HasSubstructMatch(ribose_pattern) and not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No compatible sugar ring (ribose or deoxyribose) found"

    # Check for nitrogenous bases or any common nucleotide bases
    base_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),
        Chem.MolFromSmarts("n1cnc2c1ncnc2N"),  # Adenine
        Chem.MolFromSmarts("Nc1ncnc2ncnn12"),  # Guanine
        Chem.MolFromSmarts("Nc1ncnc2ncnc2n1"),  # Hypoxanthine (for example)
        Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1=O"),  # Uracil
        Chem.MolFromSmarts("C=1N=C(NC=N1)N"),  # Cytosine
    ]
    
    base_match_found = any(mol.HasSubstructMatch(base_pattern) for base_pattern in base_patterns)
    if not base_match_found:
        return False, "No nitrogenous base found in the structure"

    return True, "Contains nucleoside with phosphate, sugar, and base attached appropriately"