"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:50362 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(O)=O)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for hydroxy group (-OH) attached to C2 (alpha position)
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[CX4H2]([CX3](=O)[OX1H0-])[OX2H1]")
    if mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return True, "Molecule contains a hydroxy group in the alpha (2) position"

    # Look for hydroxy group (-OH) attached to C2 (beta position)
    beta_hydroxy_pattern = Chem.MolFromSmarts("[CX4H2]([CX4][CX3](=O)[OX1H0-])[OX2H1]")
    if mol.HasSubstructMatch(beta_hydroxy_pattern):
        return True, "Molecule contains a hydroxy group in the beta (2) position"

    # Check for long carbon chain (fatty acid)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found (not a fatty acid)"

    return False, "No 2-hydroxy group found on a fatty acid"