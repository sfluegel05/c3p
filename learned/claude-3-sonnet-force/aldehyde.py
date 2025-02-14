"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:26718 aldehyde
An aldehyde is a compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    # Check if any carbonyl group meets the aldehyde criteria
    for match in carbonyl_matches:
        carbonyl_atom = mol.GetAtomWithIdx(match[0])
        if is_valid_aldehyde_group(mol, carbonyl_atom):
            return True, "Contains an aldehyde group (C(=O)H)"

    return False, "No valid aldehyde group found"

def is_valid_aldehyde_group(mol, carbonyl_atom):
    """
    Checks if a given carbonyl group meets the criteria for an aldehyde.

    Args:
        mol (Mol): RDKit molecule object
        carbonyl_atom (Atom): Atom object representing the carbonyl carbon

    Returns:
        bool: True if the carbonyl group is a valid aldehyde group, False otherwise
    """
    # Check if the carbonyl carbon has exactly one hydrogen neighbor
    hydrogen_count = sum(1 for neighbor in carbonyl_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
    if hydrogen_count != 1:
        return False

    # Check if the carbonyl carbon has exactly one other neighbor
    other_neighbor_count = sum(1 for neighbor in carbonyl_atom.GetNeighbors() if neighbor.GetAtomicNum() != 1 and neighbor.GetAtomicNum() != 8)
    if other_neighbor_count != 1:
        return False

    return True