"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings
    joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific corrin macrocycle pattern
    corrin_pattern = Chem.MolFromSmarts("[Co]1234N5C6=C(C)C7=[N+]1C(=CC8=[N+]2C(=C(C)C9=[N+]3[C@@](C)([C@@]5([H])[C@H](CC(N)=O)[C@@]6(C)CCC(=O)NC[C@@H](C)OP([O-])(=O)O[C@H]3[C@@H](O)[C@H](O[C@@H]3CO)n3c[n+]4c4cc(C)c(C)cc34)[C@@](C)(CC(N)=O)[C@@H]9CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]8CCC(N)=O)C(C)(C)[C@@H]7CCC(N)=O")
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "No corrin macrocycle found"

    # Check for the presence of four pyrrole-like rings (reduced or partly reduced)
    pyrrole_pattern = Chem.MolFromSmarts("[nH]1cccc1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4:
        return False, f"Found {len(pyrrole_matches)} pyrrole-like rings, need at least 4"

    # Check for the presence of three =C- groups within the macrocycle
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 3:
        return False, f"Found {len(double_bond_matches)} =C- groups, need at least 3"

    # Check for the presence of a direct carbon-carbon bond linking alpha positions within the macrocycle
    carbon_carbon_pattern = Chem.MolFromSmarts("[C]-[C]")
    carbon_carbon_matches = mol.GetSubstructMatches(carbon_carbon_pattern)
    if len(carbon_carbon_matches) < 1:
        return False, "No direct carbon-carbon bond found"

    return True, "Contains a corrin macrocycle with four pyrrole-like rings, three =C- groups, and one direct carbon-carbon bond"