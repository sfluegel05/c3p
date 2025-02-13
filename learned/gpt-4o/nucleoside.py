"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside is an N-glycosyl compound with a nucleobase attached to either a ribose or deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for nucleobases
    nucleobases_patterns = [
        Chem.MolFromSmarts("c1ncnc2c1ncnc2"),  # general purine pattern (matches adenine, guanine)
        Chem.MolFromSmarts("c1ncnc2c1nc[nH]c2"),  # pyrimidine pattern (matches cytosine, uracil, thymine)
    ]
    
    # Define SMARTS pattern for ribose (deoxyribose differs by the absence of one hydroxyl)
    ribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](C(O)CO)O1")
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](CO)C1")

    # Check for presence of a nucleobase
    has_nucleobase = any(mol.HasSubstructMatch(pat) for pat in nucleobases_patterns)
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Check for ribose or deoxyribose
    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)

    if not (has_ribose or has_deoxyribose):
        return False, "No sugar (ribose or deoxyribose) found"

    # Check for absence of phosphate group (no phosphates in nucleosides)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group(s) found, likely a nucleotide"
    
    return True, "Contains a nucleobase linked to a ribose or deoxyribose sugar"

# Test the function with examples
examples = [
    "O1[C@H](N2C=C(C(=O)NC2=O)C(O)=O)C[C@H](O)[C@H]1CO",  # Valid nucleoside
    "ClC1=NC(N)=C2N=CN(C2=N1)[C@@H]3O[C@H](COS(=O)(=O)N)[C@H]([C@H]3O)O",  # Nucleotide (should not classify as nucleoside)
]

for smi in examples:
    result, reason = is_nucleoside(smi)
    print(f"SMILES: {smi} -> Is nucleoside? {result}. Reason: {reason}")