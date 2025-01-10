"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid is characterized by a carboxyl group at position 1
    and a hydroxyl group at the omega (last) position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxyl group ('C(=O)O') at the beginning
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches or carboxyl_matches[0][0] != 0:
        return False, "No terminal carboxyl group found"

    # Identify potential terminal hydroxyl (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[C]O")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    # Determine the terminal hydroxyl, ensure it is at the end of the longest chain
    longest_chain_length = rdMolDescriptors.CalcNumAtoms(mol, onlyHeavy=True)
    terminal_hydroxyl = False
    for match in hydroxyl_matches:
        if match[0] == longest_chain_length - 1:  # Check if it's at the end
            terminal_hydroxyl = True
            break

    if not terminal_hydroxyl:
        return False, "No terminal omega-hydroxyl group found"

    # Ensure linearity: check if the hydrocarbon chain is mainly linear
    # This assumes the molecule must have a primary chain longer than a minimal size
    if longest_chain_length < 8:  # A usual minimum for definition
        return False, "Chain length too short for omega fatty acid"

    return True, "Matches omega-hydroxy fatty acid structure"