"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbon atoms in potential aliphatic chain
    chain_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Check that there is a significant aliphatic chain length
    if len(chain_atoms) < 4:
        return False, "Carbon chain length too short to be a fatty acid"

    # Check for aromaticity
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic structures"

    # Allow minimal cyclic structures only if the carboxylic group is connected to a larger chain
    ring_info = mol.GetRingInfo()
    cyclic = any(len(path) > 10 for path in ring_info.AtomRings())  # Allow small functionalities
    if cyclic:
        return False, "Significant cyclic structures are present"

    # Check if any nitrogen atoms indicate the peptide or amino linkage
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Contains amine or peptide linkages"

    # Count any disqualifying functional groups
    disallowed_functional_patterns = [
        Chem.MolFromSmarts("[S]"),  # Sulfur atoms
        Chem.MolFromSmarts("[P]"),  # Phosphoric atoms
        Chem.MolFromSmarts("[Si]")  # Silicon for large inorganics
    ]
    
    if any(mol.HasSubstructMatch(pattern) for pattern in disallowed_functional_patterns):
        return False, "Contains disallowed atom types."

    return True, "Contains aliphatic carbon chain with a carboxylic acid group, fits fatty acid profile"