"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine where the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) indicating if it's an aralkylamine and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find an amine (N with attached hydrogens possibly attached to C, excluding amides)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No primary, secondary, or tertiary amine group found"

    # Identify aromatic atoms
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    if not aromatic_atoms:
        return False, "No aromatic group found"

    # Look for connections: Amine to alkyl chain to aromatic group
    # A simple aromatic group connected through a flexible alkyl chain to an amine
    aralkylamine_pattern = Chem.MolFromSmarts("[NX3][C,c][C,c]*[a]")  # Broad N-C-C-aromatic pattern
    if mol.HasSubstructMatch(aralkylamine_pattern):
        return True, "Molecule is an aralkylamine"

    return False, "No alkyl group substituted by aromatic group found attached to amine."