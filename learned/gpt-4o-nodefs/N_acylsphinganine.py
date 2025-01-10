"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.

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

    # Check for refined sphinganine backbone with chirality and specific configuration
    sphinganine_backbone = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)[C@@H](CCCCCCCCCCCC)")

    if not mol.HasSubstructMatch(sphinganine_backbone):
        return False, "No refined sphinganine backbone found"

    # Check for amide linkage with long aliphatic chain
    acyl_chain_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)CCCCCCCCCC")
    
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No N-acyl linkage with adequate chain found"

    # Optional: Check for possible headgroup moiety if needed for specificity
    # Headgroup pattern example, may include specific sugars or other moieties
    # headgroup_pattern = Chem.MolFromSmarts("[C@H](O)CO")

    # if not mol.HasSubstructMatch(headgroup_pattern):
    #     return False, "No expected headgroup moiety found"

    # If matches all characteristics of an N-acylsphinganine
    return True, "Contains refined features consistent with N-acylsphinganine"