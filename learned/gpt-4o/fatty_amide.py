"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define amide group pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    
    # Check for presence of amide group
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"
    
    # Detect carbon chain length
    c_count = sum(1 for atom in mol.GetAtoms()
                  if atom.GetAtomicNum() == 6 and atom.GetDegree() > 1)  # Avoid terminal CH3 groups
    
    # Threshold for fatty acids (typically C12 or longer)
    if c_count < 12:
        return False, "Carbon chain too short to be a fatty acid-derived amide"
    
    return True, "Contains an amide group and a long carbon chain characteristic of fatty amides"