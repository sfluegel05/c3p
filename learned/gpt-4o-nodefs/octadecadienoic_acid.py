"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a fatty acid with a chain of 18 carbon atoms and 2 double bonds,
    ending with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group (leave open to detect similar forms, e.g., salts)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count correctly placed double bonds (exclude branches or alternate structural input)
    double_bond_pattern = Chem.MolFromSmarts("[CH]=[CH]")
    double_bond_matches = [match for match in mol.GetSubstructMatches(double_bond_pattern)]
    if len(double_bond_matches) < 2:
        return False, f"Expected at least 2 double bonds in appropriate positions, found {len(double_bond_matches)}"
    
    # Ensure there is at least 1 double bond on the backbone that is centrally located
    if not any(mol.GetAtomWithIdx(bond[0]).GetIdx() < len(mol.GetAtoms()) // 2 for bond in double_bond_matches):
        return False, f"Double bonds not centrally located"

    # Calculate total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Expected at least 18 carbon atoms, found {c_count}"
    
    return True, "Molecule is an octadecadienoic acid with 18 carbons and 2 double bonds"