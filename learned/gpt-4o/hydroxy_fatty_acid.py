"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is defined as a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the presence of a carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for the presence of at least one hydroxy group (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy (OH) group found"
    
    # Check for sufficiently long carbon chains (aliphatic)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Assuming a simple minimum chain length (can be adjusted)
        return False, "Too few carbon atoms for a fatty acid"

    # Check if chain is mostly linear (e.g., limited branching/cycles significant in fatty acids)
    # Fatty acids are primarily linear; however, they can have some branches or minor rings
    if not mol.GetRingInfo().IsAtomInRingOfSize(6): 
        return True, "Contains a carboxylic acid group and one or more hydroxy groups; valid fatty acid configuration"
    else:
        return False, "Presence of cyclic structures inconsistent with a typical fatty acid"

    return True, "Contains a carboxylic acid group and one or more hydroxy groups on a sufficiently long carbon chain"