"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Identify if there is a terminal hydroxy group at omega position
    omega_hydroxy_pattern = Chem.MolFromSmarts("OCC(=O)O")
    if not mol.HasSubstructMatch(omega_hydroxy_pattern):
        return False, "No omega-terminal hydroxy group found"

    # Check for a carboxylic acid group in the molecule
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Verify the length of the carbon chain (minimum 7 for a short fatty acid)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 7:
        return False, f"Carbon chain too short, found {len(carbon_atoms)} carbons"

    return True, "Contains omega-terminal hydroxy and carboxylic acid groups with sufficient chain length"