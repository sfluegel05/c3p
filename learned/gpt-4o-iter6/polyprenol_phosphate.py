"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate consists of a polyisoprene chain with a terminal
    allylic hydroxy group condensed with a phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a polyprenol phosphate, False otherwise.
        str: Reason for classification or failure.
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Isoprene pattern with possible cis/trans configurations
    isoprene_pattern = Chem.MolFromSmarts("C(=C)CC")
    
    # Polyisoprenoid chain pattern: series of isoprene units
    polyisoprenoid_pattern = Chem.MolFromSmarts("C(=C)(CC){5,}") # Require at least 5 isoprene units
    
    # Phosphate group patterns
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)[O-]")
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])OP(=O)(O)[O-]")

    # Check the presence of long polyisoprenoid chain
    if not mol.HasSubstructMatch(polyisoprenoid_pattern):
        return False, "No sufficient polyisoprenoid chain detected"

    # Check phosphate attachment
    phosphate_matches = mol.HasSubstructMatch(phosphate_pattern)
    diphosphate_matches = mol.HasSubstructMatch(diphosphate_pattern)

    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate or diphosphate group found"

    # Ensure that phosphate/diphosphate attaches to an allylic carbon
    allylic_phosphate_pattern = Chem.MolFromSmarts("C[O]P(=O)(O)O")
    if not mol.HasSubstructMatch(allylic_phosphate_pattern):
        return False, "Phosphate group not correctly attached to an allylic position"

    return True, "Valid polyprenol phosphate"