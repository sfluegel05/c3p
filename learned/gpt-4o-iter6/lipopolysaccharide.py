"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain a specific trisaccharide unit, oligosaccharide side chains, and long fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific and general patterns needed
    # Hypothetical SMARTS for simplicity:
    trisaccharide_pattern_1 = Chem.MolFromSmarts("[C@H]([C@H](O)[C@H](O)[C@H]O)[C@H](O)C(O)=O") # Simplified unit placeholder
    trisaccharide_pattern_2 = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H](O)[C@H](O1)O)C(O)=O") # Another simplified pattern
    
    if not (mol.HasSubstructMatch(trisaccharide_pattern_1) or mol.HasSubstructMatch(trisaccharide_pattern_2)):
        return False, "No specific trisaccharide units found"

    # Consider variability in the fatty acid motifs - example hydroxyl derivatives
    long_fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCC(CC)CCCC(=O)O") # Example adjusted for specificity
    if not mol.HasSubstructMatch(long_fatty_acid_pattern):
        return False, "No long fatty acid chains detected"

    # Broadening oligosaccharide scope
    oligosaccharide_pattern = Chem.MolFromSmarts("[C@H](CO)O[C@H][C@H](O)O") # Example broader representative
    if not mol.HasSubstructMatch(oligosaccharide_pattern):
        return False, "No sufficient oligosaccharide side chains detected"

    # Complex molecules characteristic patterns: rings, oxygen content, and bonds
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_rings < 4 or oxygen_count < 20:
        return False, "Insufficient complexity or polysaccharide indication"

    return True, "Matches key structural features of lipopolysaccharides: complex units, oligosaccharides, and fatty acids"