"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the deprotonated form of a fatty acid, characterized by a carboxylate group (-COO-) at the end of a long hydrocarbon chain.

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
    
    # Check for a terminal carboxylate group connected to a long, non-branching carbon chain (6 to 22 carbons)
    # We use a more specific SMARTS pattern to ensure the carboxylate is at the end of the chain and
    # that we have a long, unbranched chain connected to the carboxylate carbon.
    # We check for chains ranging from 6 to 22 carbons.
    for chain_length in range(6, 23):
      chain_smarts = f"[C](=[O])[O-]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-" + "-[CH2]"*(chain_length-6) + "[CH3,CH2]"
      pattern = Chem.MolFromSmarts(chain_smarts)
      if mol.HasSubstructMatch(pattern):
        return True, f"Contains terminal carboxylate group and a long, unbranched carbon chain of length >= {chain_length}, which is consistent with a fatty acid anion."

    return False, "Does not contain a carboxylate group at the end of a long carbon chain with at least 6 carbons."