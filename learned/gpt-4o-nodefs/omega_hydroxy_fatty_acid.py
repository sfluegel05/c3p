"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    Omega-hydroxy fatty acids have a hydroxyl group at the terminal (omega) position and a carboxylic acid group.

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

    # Identify omega-terminal hydroxy group
    terminal_hydroxy_pattern = Chem.MolFromSmarts("[CH2]CO")
    if not mol.HasSubstructMatch(terminal_hydroxy_pattern):
        return False, "No omega-terminal hydroxy group found"

    # Check for presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Ensure sufficient chain length for fatty acid
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 7:
        return False, f"Carbon chain too short, found {len(carbon_atoms)} carbons"

    return True, "Contains omega-terminal hydroxy and carboxylic acid groups with sufficient chain length"