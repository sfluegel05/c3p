"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:35526 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as an amide of a sulfonic acid: RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonamide pattern S(=O)(=O)N
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    matches = mol.GetSubstructMatches(sulfonamide_pattern)
    
    if not matches:
        return False, "Does not contain sulfonamide group S(=O)(=O)N"

    # Check for common substituents on sulfur and nitrogen
    s_atoms = [mol.GetAtomWithIdx(match[0]).GetNeighbors() for match in matches]
    n_atoms = [mol.GetAtomWithIdx(match[1]).GetNeighbors() for match in matches]
    
    valid_s_substituents = [atom for atom in s_atoms if any(neighbor.GetAtomicNum() == 6 for neighbor in atom)]
    valid_n_substituents = [atom for atom in n_atoms if any(neighbor.GetAtomicNum() == 6 for neighbor in atom)]

    if valid_s_substituents and valid_n_substituents:
        return True, "Contains sulfonamide group S(=O)(=O)N with valid substituents"
    else:
        return False, "Contains S(=O)(=O)N group but lacks valid substituents"