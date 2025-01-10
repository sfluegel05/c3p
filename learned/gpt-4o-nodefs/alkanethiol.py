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
    thiol_pattern = Chem.MolFromSmarts("[#6][SX2H]")  # Carbon connected to SH
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No correct thiol groups found"
    
    # Check for presence of aromatic rings
    if mol.GetAromaticAtoms():
        return False, "Contains aromatic rings"

    # Allowable elements in alkanethiol
    allowable_atoms = {'C', 'H', 'S', 'O', 'N'}

    # Count the atoms and ensure only allowable elements are present
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowable_atoms:
            return False, f"Contains non-alkanethiol characteristic atoms: {atom.GetSymbol()}"
    
    # Check the number of potential oxygens which should not lead to ester formation
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if o_count > 4:
        return False, "Too many oxygen atoms suggesting non-alkanethiol complexity"

    return True, "Contains thiol group with alkane or simple aliphatic association"