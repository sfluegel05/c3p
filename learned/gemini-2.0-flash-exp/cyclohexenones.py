"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a 6 membered ring. Substituents are allowed
    ring_pattern = Chem.MolFromSmarts("[C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3]")
    if not mol.HasSubstructMatch(ring_pattern):
      return False, "Molecule does not contain a six membered ring"
        
    # Define SMARTS pattern for carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("[C;R1]=[O]")
    
    # Define SMARTS pattern for double bond
    double_bond_pattern = Chem.MolFromSmarts("[C;R1]=[C;R1]")
    
    # Check if molecule contains both a carbonyl and a double bond in the same ring
    if mol.HasSubstructMatch(carbonyl_pattern) and mol.HasSubstructMatch(double_bond_pattern):
      return True, "Molecule is a cyclohexenone"
    
    return False, "Molecule does not contain both a carbonyl and a double bond in the same ring"