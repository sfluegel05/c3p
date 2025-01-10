"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    Tetraterpenoids typically derive from a C40 backbone with modifications, and feature a polyene structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check number of carbon atoms, aiming for a backbone capable of modification (typically ~40)
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 35 or num_carbons > 50:  # Broadening range for flexibility
        return False, f"Incorrect number of carbons for tetraterpenoid: {num_carbons} carbons"

    # Check presence of multiple double bonds characteristic of a polyene
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() >= 2)
    if num_double_bonds < 10:  # Expecting many double bonds for polyene system
        return False, "Insufficient double bonds for a tetraterpenoid structure"
    
    # Check for isoprene-like repeating units (~8 isoprene units expected)
    isoprene_unit = Chem.MolFromSmarts("C=C-C-C=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit)
    if len(isoprene_matches) < 5:  # Lowering threshold to account for modifications/rearrangements
        return False, "Not enough isoprene-like units to suggest tetraterpenoid origin"
    
    # If all checks pass
    return True, "Structure likely contains the characteristic features of a tetraterpenoid"