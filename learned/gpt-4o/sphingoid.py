"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids generally have a long carbon chain, an amino group, and hydroxyl groups,
    as well as potential additional features like double bonds or certain stereochemistry.
    
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

    # Check for a long aliphatic carbon chain (14-30+ carbons)
    long_chain_pattern = Chem.MolFromSmarts("[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficient long carbon chain recognized"

    # Check for the presence of a 2-amino-1,3-diol configuration
    sphingoid_core_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](N)CO")
    if not mol.HasSubstructMatch(sphingoid_core_pattern):
        return False, "No recognizable sphingoid core substructure"

    # Check for additional hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Less than two hydroxyl groups found"

    # Evaluate the total number of carbon atoms
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 14:
        return False, f"Carbon chain length {total_carbons} is insufficient for a sphingoid"

    # Check for common side features like double bonds or extra hydroxyls
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains additional features like double bonds, consistent with sphingoid derivatives"

    return True, "Contains typical structure of a sphingoid (long chain, amino, hydroxyl groups)"