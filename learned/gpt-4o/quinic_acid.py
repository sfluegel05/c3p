"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or one of its derivatives based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid characterized by a cyclohexane backbone with an
    attached carboxylic acid group and several hydroxyl groups, possibly with further esterification.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as quinic acid or a derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Corrected pattern for hydroxylated cyclohexane with carboxylic acid
    quinic_acid_base_pattern = Chem.MolFromSmarts("[C@@H]1[C@@H](O)[C@@H](O)[C@@H](O)C[C@@H]1C(=O)O")
    
    # Check for cyclohexane ring with required stereochemistry and functional groups
    if not mol.HasSubstructMatch(quinic_acid_base_pattern):
        return False, "No cyclohexane with correct stereochemistry and hydroxyl groups found plus carboxylic acid group"
    
    # Optional: Check for additional ester or ether linkages
    ester_linkage_pattern = Chem.MolFromSmarts("[C](=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1")
    if mol.HasSubstructMatch(ester_linkage_pattern):
        return True, "Quinic acid derivative with ester linkage detected"
    
    # Additional indication of complex esterification/ether linkages (consider aromatic systems)
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for smarts in [
        "[c][C](=O)O[C@@H]1C[C@@](O)(C[C@@H](O)[C@@H]1O)C(=O)O",
        "[c][C](=O)O[C@@H]1C[C@@H](O)[C@@](O)(C[C@@H]1O)C(=O)O"
    ]):
        return True, "Complex quinic acid ester derivative"

    return True, "Basic quinic acid structure"