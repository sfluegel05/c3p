"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to accurately count hydrogen atoms
    mol = Chem.AddHs(mol)

    # Initialize flag
    is_alkanethiol = False

    # Iterate over atoms to find sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Count hydrogens attached to sulfur
            num_hydrogens = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
            # Check if sulfur has exactly one hydrogen attached
            if num_hydrogens == 1:
                # Check for carbon neighbors
                has_carbon_neighbor = any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors())
                # Check for other heteroatoms connected to sulfur
                has_other_heteroatoms = any(
                    neighbor.GetAtomicNum() not in [1, 6, 16] for neighbor in atom.GetNeighbors()
                )
                if has_carbon_neighbor and not has_other_heteroatoms:
                    is_alkanethiol = True
                    break  # Found the required group

    if is_alkanethiol:
        return True, "Contains a sulfanyl group (-SH) attached to an alkyl group"
    else:
        return False, "Does not contain an -SH group attached to an alkyl group"