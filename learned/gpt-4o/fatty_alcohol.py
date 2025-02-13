"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"
    
    # Check for potential aromatics but don't dismiss outright
    aromatic_count = sum(atom.GetIsAromatic() for atom in mol.GetAtoms())
    total_atoms = mol.GetNumHeavyAtoms()
    if aromatic_count > 0.3 * total_atoms:  # Allow for up to 30% aromatic
        return False, "Molecule is too aromatic"

    # Check for existence of suitable carbon chains (including unsaturation)
    carbon_chain_pattern = Chem.MolFromSmarts("CCCC")  # Simple flexible carbon chain, possibly to be expanded
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No suitable aliphatic carbon chain found"

    # Count the number of carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Fatty alcohols need at least 3 carbons and 1 hydroxyl group
    if c_count < 3 or o_count < 1:
        return False, "Too few carbon atoms or hydroxyl groups for a fatty alcohol"

    # Ensure it is majorly aliphatic (considering branching and lack of complex aromatic rings)
    if c_count / total_atoms < 0.75:  # Modified threshold to consider branching
        return False, f"Structure is not primarily aliphatic with {c_count} out of {total_atoms} heavy atoms being carbon"

    return True, f"Contains suitable alcohol group with {c_count} carbon atoms forming an aliphatic chain"