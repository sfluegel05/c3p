"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group (-C([O-])=O)
    and a hydrocarbon chain, which may include double bonds or functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Require a carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group (-C([O-])=O) found"

    # Check for a chain with at least 6 continuous carbon atoms (allowing for branching, double bonds, or heteroatoms)
    # Simplify the pattern to detect long carbon chains, given the variety of structures
    carbon_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C")  # Matches a chain of any carbon connectivity
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient carbon chain length for fatty acid anion"

    return True, "Contains a deprotonated carboxylate group with a sufficiently long carbon chain"