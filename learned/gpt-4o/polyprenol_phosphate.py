"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate includes a polyprenol chain attached via a phosphate ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphate or diphosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("[OX1][PX4](=[OX1])(O)")
    diphosphate_pattern = Chem.MolFromSmarts("[OX1][PX4](=[OX1])([OX1])O[PX4](=[OX1])")
    if not (mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(diphosphate_pattern)):
        return False, "No phosphate or diphosphate group found"
    
    # Check for ester linkage to phosphate group: O-P(=O)(O)...
    ester_linkage_smarts = "O[PX4](=O)(O)"
    ester_linkage = Chem.MolFromSmarts(ester_linkage_smarts)
    if not mol.HasSubstructMatch(ester_linkage):
        return False, "Ester linkage to phosphate group not found"

    # Look for repeating isoprene units (better pattern for polyprenol) in linear fashion
    # We are enhancing to identify C=C connectivity more precisely.
    isoprene_pattern = Chem.MolFromSmarts("C=CC(CC=C)")
    isoprene_count = len(mol.GetSubstructMatches(isoprene_pattern))
    
    # It's a polyprenol phosphate if there are 3 or more.
    if isoprene_count < 3:
        return False, f"Insufficient isoprene units: {isoprene_count}"
    
    return True, "Contains polyprenol chain with phosphate ester bond"