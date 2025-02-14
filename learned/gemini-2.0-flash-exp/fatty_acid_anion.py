"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the deprotonated form of a fatty acid, characterized by a carboxylate group (-COO-) and a hydrocarbon chain.

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
    
    # Check for carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
       return False, "No carboxylate group found"

    # Check for at least three carbons directly attached to the carboxylate carbon
    alkyl_pattern = Chem.MolFromSmarts("C(=O)([O-])[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(alkyl_pattern):
      return False, "Carboxylate group is not attached to sufficient carbon chain"
        
    # Check for a carbon chain of at least 3 carbons
    # this could be improved by checking a specific length but for now this is good
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
      return False, "No sufficient carbon chain detected"

    
    return True, "Contains carboxylate group and attached carbon chain (at least 3 carbons), which is consistent with a fatty acid anion."