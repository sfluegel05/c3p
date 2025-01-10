"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a long hydrocarbon chain and a deprotonated carboxylic acid group (-C([O-])=O).

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

    # SMARTS pattern for carboxylate group -C([O-])=O
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group (-C([O-])=O) found"
    
    # Check for a sufficiently long carbon chain
    # Arbitrarily define "long" as having at least 7 carbon atoms apart from the carboxylate
    chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[O-]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Too few carbon atoms for fatty acid chain"

    # If both patterns are found, we classify it as a fatty acid anion
    return True, "Contains a deprotonated carboxylate group with a long hydrocarbon chain"