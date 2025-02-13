"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group (-C([O-])=O)
    and a hydrocarbon chain, which may include double bonds, but typically form straight or slightly branched chains.

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

    # Define carboxylate SMARTS pattern
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if carboxylate_pattern is None:
        return (None, "Error creating carboxylate pattern")

    # Detect the presence of a carboxylate group
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group (-C([O-])=O) found"

    # Define a pattern for a simple hydrocarbon chain
    # Allow straight or lightly branched alkyl chains
    chain_pattern = Chem.MolFromSmarts("[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]")
    if chain_pattern is None:
        return (None, "Error creating chain pattern")

    # Ensure it contains a sufficiently long chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Hydrocarbon chain does not match typical fatty acid anion structure"

    return True, "Contains a deprotonated carboxylate group with a suitable hydrocarbon chain characteristic of fatty acid anions"