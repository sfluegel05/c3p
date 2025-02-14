"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has a carboxylic acid group, a long carbon chain and one or more alkyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a chain of at least 3 carbons
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]([CX4,CX3])[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No 3-carbon chain found"
    
    # Check for a branching point
    branch_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6: # Check only carbon atoms
            neighbors = [neighbor for neighbor in atom.GetNeighbors()]
            carbon_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() == 6]
            non_h_neighbors = [neighbor for neighbor in neighbors if neighbor.GetAtomicNum() != 1]
            
            if len(non_h_neighbors) > 2 and len(carbon_neighbors) >= 2:
                branch_found = True
                break

    if not branch_found:
      return False, "No appropriate branching point found"

    # Count carbons to ensure minimum length.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
      return False, "Molecule is too short"
      
    return True, "Contains a carboxylic acid group, a long chain and one or more branching points"