"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    Tetraterpenoids typically have a C40 backbone with various possible modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check number of carbon atoms, typically 40 for tetraterpenoids
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 35 or num_carbons > 42:  # allowing some flexibility
        return False, f"Incorrect number of carbons for tetraterpenoid: {num_carbons} carbons"

    # Check presence of unsaturated bonds, a characteristic of terpenes
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() >= 2)
    if num_double_bonds < 8:  # assuming unsaturation presence
        return False, "Insufficient double bonds for a tetraterpenoid structure"
    
    # Detect common modifications: alcohol or carbonyl groups
    alcohol_pattern = Chem.MolFromSmarts('[OX2H]')  # Hydroxyl group
    carbonyl_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')  # Carbonyl group (>C=O)

    has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    if not (has_alcohol or has_carbonyl):
        return False, "Lacks typical tetraterpenoid modification groups (alcohol or carbonyl)"
    
    # If all checks pass
    return True, "Structure likely contains the characteristic features of a tetraterpenoid"