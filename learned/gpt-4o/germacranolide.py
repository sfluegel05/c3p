"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is defined as a sesquiterpene lactone with a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for germacrane backbone (approximated, focusing on key features typically found)
    germacrane_pattern = Chem.MolFromSmarts("C1=CCCC[C@H]2C[C@H](C)C=C/C(=O)O[C@@H]12")
    
    # Define SMARTS pattern for lactone group (-C(=O)O-)
    lactone_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for germacrane-like backbone pattern
    if not mol.HasSubstructMatch(germacrane_pattern):
        return False, "No germacrane-like structure found"
    
    # Ensure it contains a lactone group
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"
    
    # Check specific connectivity traits for sesquiterpene characteristics
    double_bonds_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bonds_pattern)
    if len(double_bond_matches) < 2:  # at least two implicit or explicit double bonds
        return False, "Insufficient double bonds for germacranolide"

    # Ensure the presence of oxygen functionalities (like hydroxyls commonly present)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Too few oxygens for typical germacranolide functionalization"

    # Further checks such as stereochemistry validations could be performed if needed
    # Additional check on rotatable bonds found in typical germacranolides
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Too few rotatable bonds for typical sesquiterpene complexity"

    return True, "Contains germacrane skeleton and key functional groups typical of germacranolides"