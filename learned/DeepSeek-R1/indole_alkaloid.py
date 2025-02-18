"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:24828 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must contain an indole skeleton (bicyclic structure with 
    benzene fused to pyrrole ring) and have at least one additional nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define indole core pattern (benzene fused to pyrrole)
    indole_pattern = MolFromSmarts("n1ccc2ccccc12")
    
    # Check for indole skeleton
    indole_matches = mol.GetSubstructMatches(indole_pattern)
    if not indole_matches:
        return False, "No indole skeleton found"

    # Collect all atoms in indole substructures
    indole_atoms = set()
    for match in indole_matches:
        indole_atoms.update(match)

    # Check for at least one nitrogen outside the indole skeleton
    has_external_nitrogen = any(
        atom.GetAtomicNum() == 7 and atom.GetIdx() not in indole_atoms
        for atom in mol.GetAtoms()
    )

    if not has_external_nitrogen:
        return False, "No additional nitrogen outside indole skeleton"

    return True, "Contains indole skeleton with additional nitrogen characteristic of alkaloids"