"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    Catecholamines typically have a catechol moiety (benzene with hydroxyl groups)
    and an amine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for catechol moiety (benzene with hydroxyl groups, allowing flexibility)
    catechol_moiety = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1 |c1cc(O)c(O)cc1|")  # Allow flexibility in hydroxyl positioning
    catechol_matches = mol.GetSubstructMatches(catechol_moiety)
    if not catechol_matches:
        return False, "No catechol moiety found"
    
    # Look for an amine group (allowing tertiary and quaternary amines)
    amine_patterns = [
        Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]"),  # Terminal nitrogen patterns (primary, secondary, tertiary amines)
        Chem.MolFromSmarts("[NX4+]")  # Quaternary amine (charged nitrogen)
    ]

    amine_found = any(mol.HasSubstructMatch(pattern) for pattern in amine_patterns)
    if not amine_found:
        return False, "No amine group found"
    
    return True, "Contains catechol moiety and amine group indicating a catecholamine"

# Example usage
print(is_catecholamine("CC(C)NC[C@H](O)c1ccc(O)c(O)c1"))  # L-isoprenaline
print(is_catecholamine("COC1=CC=CC=C1"))  # Anisole (should not be catecholamine)