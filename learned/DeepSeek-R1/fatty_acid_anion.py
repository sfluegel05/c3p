"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI: fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, characterized by a carboxylate group (-COO-).

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

    # Check for carboxylate group (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts('[O-]C(=O)')
    if mol.HasSubstructMatch(carboxylate_pattern):
        return True, "Contains carboxylate group (deprotonated carboxylic acid)"
    else:
        return False, "No carboxylate group found"