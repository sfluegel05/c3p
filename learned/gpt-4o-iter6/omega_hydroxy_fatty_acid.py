"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl group at the last position along a linear chain.

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

    # Find the longest chain in the structure
    longest_chain = Chem.rdmolops.GetLongestPath(mol)
    if not longest_chain or len(longest_chain) < 8:
        return False, "The main hydrocarbon chain is too short (<8 carbons)."

    # Check if there's a carboxyl group (COOH) at position 1
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches or carboxyl_matches[0][0] != longest_chain[0]:
        return False, "No terminal carboxyl group found at the start of the chain."

    # Check if there's a hydroxyl group (OH) at the last position
    hydroxyl_pattern = Chem.MolFromSmarts("CO")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    terminal_hydroxyl = any(match[1] == longest_chain[-1] for match in hydroxyl_matches)

    if not terminal_hydroxyl:
        return False, "No terminal omega-hydroxyl group found at the end of the chain."

    return True, "Molecule is an omega-hydroxy fatty acid: has terminal carboxyl and omega-hydroxyl groups on a linear chain."