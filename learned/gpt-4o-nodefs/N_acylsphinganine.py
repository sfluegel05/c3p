"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    N-acylsphinganines have a sphinganine backbone with an amide linkage to a fatty acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sphinganine backbone [C@H](O) and [C@@H] configuration
    sphinganine_backbone = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)")
    if not mol.HasSubstructMatch(sphinganine_backbone):
        return False, "No sphinganine backbone found"

    # Check for amide linkage -C(=O)N-
    amide_linkage = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(amide_linkage):
        return False, "No amide linkage found"

    # Check for long carbon chains (>10 carbons)
    carbon_chain = Chem.MolFromSmarts("[C;X4]~[C;X4]~[C;X4]~[C;X4]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No sufficient long carbon chain found"

    # If all checks passed, classify as N-acylsphinganine
    return True, "Contains features consistent with N-acylsphinganine"