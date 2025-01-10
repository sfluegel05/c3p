"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find potential oligosaccharide patterns (a general sugar pattern)
    sugar_pattern = Chem.MolFromSmarts("[Cx4]1[Ox1][Cx4][Cx4][Ox1][Cx4]1")  # Ring structure typical in sugars
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No oligosaccharide composition found"
        
    # Find long hydrocarbon chain (fatty acid of sufficient length)
    long_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCC(C)CCCCCCCCCCC")  # Linear long chain with a carbonyl
    if not mol.HasSubstructMatch(long_fatty_acid_pattern):
        return False, "Lacks fatty acid unit"

    # Find the presence of multiple ring systems that could represent the complex structure
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3:  # Needs to have more than one ring for the trisaccharide unit
        return False, "Too few ring systems detected"

    return True, "Contains oligosaccharide structures with long fatty acids typical of lipopolysaccharides"