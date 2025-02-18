"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine with an alkyl group substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amine group (primary, secondary, or tertiary)
    amine_pattern = Chem.MolFromSmarts("N")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Look for aromatic group
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    if not aromatic_atoms:
        return False, "No aromatic group found"

    # Search for aralkylamine pattern: nitrogen linked to an aromatic via an alkyl chain
    aralkylamine_pattern = Chem.MolFromSmarts("[N]@[CX4][c]")
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "No alkyl group substituted by aromatic group attached to amine"

    return True, "Molecule is an aralkylamine"

# The function now checks for an amine group connected to an aromatic group via an alkyl chain