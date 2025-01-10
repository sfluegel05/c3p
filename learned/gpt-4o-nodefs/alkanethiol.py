"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    Alkanethiols feature one or more thiol groups (-SH) attached to an alkane or
    simple aliphatic chain excluding larger complex structures like peptides.

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
    thiol_pattern = Chem.MolFromSmarts("[SX2H]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    if len(thiol_matches) == 0:
        return False, "No thiol groups found"
    
    # Check that the molecule consists primarily of carbon with the acceptable, minimal functional atoms
    acceptable_atoms = {'C', 'H', 'S', 'O', 'N'}
    atom_counts = {atom.GetSymbol() : 0 for atom in mol.GetAtoms()}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in acceptable_atoms:
            return False, f"Contains non-alkanethiol characteristic atoms: {atom.GetSymbol()}"
        atom_counts[atom.GetSymbol()] += 1

    # Check for large counts of typically peptide/complex forming atoms: N, O 
    if atom_counts['N'] > 2 or atom_counts['O'] > 2:
        return False, "Too many functional atoms indicating non-alkanethiol complexity"
    
    # Control against complex peptide-like molecules, such as 3+ peptide bonds (C(=O)NH)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_matches) >= 3:
        return False, "Contains peptide-like characteristics"

    return True, "Contains thiol group with alkyl or simple aliphatic association"