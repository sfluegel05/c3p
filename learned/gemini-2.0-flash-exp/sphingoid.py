"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid has a long hydrocarbon chain, an amino group, and two hydroxyl or carbonyl groups, usually in vicinal position

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10 or carbon_count > 30:
       return False, f"Carbon chain length ({carbon_count}) is outside the typical range for sphingoids (10-30)"

    # Check for the presence of amino group
    amino_pattern = Chem.MolFromSmarts("[NX3,NX4]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Check for pattern with two hydroxyl groups
    diol_pattern1 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([OH1])([CX4])")
    # Check for pattern with one carbonyl and one hydroxyl group
    diol_pattern2 = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[CX4]([NX3,NX4])[CX4]")
    # Check for pattern with an additional double bond adjacent to the alcohol
    diol_pattern3 = Chem.MolFromSmarts("[CX4]=[CX3][CX4]([OH1])[CX4]([NX3,NX4])[CX4]")
    if not (mol.HasSubstructMatch(diol_pattern1) or mol.HasSubstructMatch(diol_pattern2) or mol.HasSubstructMatch(diol_pattern3)):
      return False, "Did not find the necessary hydroxyl or carbonyl groups near the amino group"

    # Check for a long chain with double bond
    chain_double_bond = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]=[CX4,CX3]~[CX4,CX3]")
    if mol.HasSubstructMatch(chain_double_bond):
        return True, "Contains long hydrocarbon chain, an amino group, and hydroxyl/carbonyl groups and double bond"

    return True, "Contains long hydrocarbon chain, an amino group, and hydroxyl/carbonyl groups"