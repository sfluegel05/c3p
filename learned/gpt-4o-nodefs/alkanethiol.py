"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    Alkanethiols feature one or more thiol groups (-SH) attached to an alkane or
    simple aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thiol group pattern (-SH)
    thiol_pattern = Chem.MolFromSmarts("[SH]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    if len(thiol_matches) == 0:
        return False, "No thiol groups found"
    
    # Ensure that the molecule is composed of acceptable alkanethiol atoms
    acceptable_atoms = {'C', 'H', 'S', 'O'}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in acceptable_atoms:
            return False, f"Contains non-alkanethiol characteristic atoms: {atom.GetSymbol()}"

    # Check for excessive oxygen atoms, often associated with non-alkanethiol groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if o_count > 2:
        return False, "Too many oxygen atoms suggesting non-alkanethiol complexity"
    
    return True, "Contains thiol group with alkane or simple aliphatic association"