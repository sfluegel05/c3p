"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine that has an alkyl group substituted by an aromatic group.

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

    # Check for aromatic groups in the molecule
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    if not aromatic_atoms:
        return False, "No aromatic group found"

    # Look for amine group (primary, secondary, or tertiary)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No primary, secondary, or tertiary amine group found"

    # Search for an aralkylamine: an amine should connect through an alkyl chain to an aromatic ring
    aralkylamine_pattern = Chem.MolFromSmarts("[N][CX4][CX4][r]")
    
    # Adjust pattern if necessary to capture more diverse structures, or add alternate patterns:
    complex_aralkylamine_pattern = Chem.MolFromSmarts("[N]([C,c])[C,c][r]")

    if not (mol.HasSubstructMatch(aralkylamine_pattern) or mol.HasSubstructMatch(complex_aralkylamine_pattern)):
        return False, "No alkyl group substituted by aromatic group found attached to amine."

    return True, "Molecule is an aralkylamine"

# This methodology attempts to utilize refined SMARTS patterns to better capture diverse aralkylamine structures.