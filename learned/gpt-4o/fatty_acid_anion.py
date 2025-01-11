"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by the presence of a carboxylate group and a long carbon chain.

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
    
    # Look for carboxylate anion pattern [O-]C(=O)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion group found"
        
    # Verify the presence of a long carbon chain (e.g., C-C-C)
    carbon_chain_pattern = Chem.MolFromSmarts("C~C~C")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found" 

    # Ensure there are more than a minimal number of carbon atoms (e.g., ≥ 10 for a minimal fatty acid)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Insufficient carbon count ({c_count}), too short to be a fatty acid anion"

    return True, "Contains carboxylate anion group with an adequate carbon chain length"