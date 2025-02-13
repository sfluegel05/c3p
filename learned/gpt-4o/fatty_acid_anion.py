"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group 
    and a significant hydrocarbon content, possibly with some functionalization.

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
    
    # Look for carboxylate group pattern [CX3](=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found (required for fatty acid anion)"
    
    # Check for a significant amount of carbon atoms for fatty acid chains
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())

    # Fatty acid anions typically have long carbon chains, we assume at least 10 carbons
    if carbon_count < 10:
        return False, f"Insufficient carbon content for a fatty acid anion (found {carbon_count} carbons)"

    # Evaluate general structural characteristics for various fatty acid anions without strict chain structure
    return True, f"Contains a carboxylate group and {carbon_count} carbon atoms common in fatty acid anions"