"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    Catecholamines have a catechol moiety (benzene with two adjacent hydroxyl groups)
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
    
    # Look for catechol moiety (benzene with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety found"
    
    # Look for an amine group
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Matches primary/secondary amines not in amides
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    return True, "Contains catechol moiety and amine group indicating a catecholamine"

# Example usage
print(is_catecholamine("CC(C)NC[C@H](O)c1ccc(O)c(O)c1"))  # L-isoprenaline
print(is_catecholamine("COC1=CC=CC=C1"))  # Anisole (should not be catecholamine)