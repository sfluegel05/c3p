"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide consists of at least two monosaccharide units joined by glycosidic linkages.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for common sugar rings (hexose and pentose sugars, different stereochemistries)
    sugar_patterns = [
        "[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O",    # Alpha-D-glucopyranose (simplified)
        "[C@@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@@H]1O",   # Beta-D-glucopyranose (simplified)
        "[C@H](O)CC(O)[C@H]1O[C@H]([C@@H](O)[C@H]1O)C", # Another common sugar pattern
    ]
    
    # Compile patterns into Mol objects
    sugar_mols = [Chem.MolFromSmarts(pat) for pat in sugar_patterns if Chem.MolFromSmarts(pat)]
    
    # Define SMARTS pattern for glycosidic linkages
    glycosidic_bond_pattern = Chem.MolFromSmarts("C1OC1")  # Simplified cyclic ether
    
    # Check for at least two sugar units
    sugar_count = 0
    for sugar_mol in sugar_mols:
        sugar_count += len(mol.GetSubstructMatches(sugar_mol))
    
    if sugar_count < 2:  # Oligosaccharides must have at least two sugar units
        return False, "Insufficient monosaccharide units"
    
    # Check for glycosidic linkages
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic linkages found"
    
    return True, "Contains sufficient monosaccharide units and glycosidic linkages indicative of an oligosaccharide"