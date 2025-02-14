"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a fatty acid containing an epoxide ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide ring pattern (more precise definition of an epoxide)
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"
    
   # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"


    # Create a SMARTS pattern for a fatty acid chain with an epoxide:
    # epoxide - chain of at least 3 C - carboxylic acid
    # The [!R] is to exclude the case where the epoxide is part of a ring
    # The [*] is a generic atom to account for single or double bonds on the chain
    fatty_acid_pattern = Chem.MolFromSmarts("[C]1[O][C]1[!R][CX4,CX3][CX4,CX3][CX4,CX3][CX4,CX3]*C(=O)O")
    
    if not mol.HasSubstructMatch(fatty_acid_pattern):
      # Check for shorter chain, at least 3 carbons
      fatty_acid_pattern = Chem.MolFromSmarts("[C]1[O][C]1[!R][CX4,CX3][CX4,CX3][CX4,CX3]*C(=O)O")
      if not mol.HasSubstructMatch(fatty_acid_pattern):
          return False, "No fatty acid chain with epoxide found"

    # Check the number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too few carbons for a fatty acid"

    return True, "Contains an epoxide ring and a fatty acid chain"