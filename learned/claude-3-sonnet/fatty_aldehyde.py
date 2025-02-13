"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has an aldehyde group at one end of a carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde group (-CH=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CH1X3](=[OX1])")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # Also check for non-explicit hydrogen aldehyde pattern
    if not aldehyde_matches:
        aldehyde_pattern = Chem.MolFromSmarts("[CX3H0](=[OX1])")
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    if len(aldehyde_matches) > 1:
        return False, "Multiple aldehyde groups found"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Carbon chain too short for fatty aldehyde"

    # Verify the aldehyde is terminal (should have only one non-H neighbor)
    aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_matches[0][0])
    non_h_neighbors = [n for n in aldehyde_carbon.GetNeighbors() 
                      if n.GetAtomicNum() != 1 and n.GetSymbol() != 'O']
    if len(non_h_neighbors) != 1:
        return False, "Aldehyde group is not terminal"

    # Look for carbon chain
    chain_pattern = Chem.MolFromSmarts("[C]~[C]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No carbon chain found"

    # Check that most atoms are carbons (allowing for some heteroatoms)
    total_atoms = mol.GetNumAtoms()
    if carbon_count < total_atoms * 0.6:
        return False, "Too many heteroatoms for fatty aldehyde"

    # Success - determine if saturated or unsaturated
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    has_double_bonds = mol.HasSubstructMatch(double_bond_pattern)
    
    if has_double_bonds:
        return True, "Unsaturated fatty aldehyde with terminal aldehyde group"
    else:
        return True, "Saturated fatty aldehyde with terminal aldehyde group"