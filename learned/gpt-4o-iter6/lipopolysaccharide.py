"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain a complex mixture of oligosaccharide structures and long fatty acid chains.

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

    # Look for trisaccharide unit pattern, simplified to sugar linkages [complex sugar pattern placeholder]
    trisaccharide_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O[C@H]1)CO")  # Minimal pattern
    if not mol.HasSubstructMatch(trisaccharide_pattern):
        return False, "No trisaccharide unit found"

    # Look for long lipid chain bound to sugar [complex lipid pattern placeholder]
    lipid_chain_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCC")  # Long hydrocarbon chain example
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No long fatty acid chain detected"
    
    # Additional complex ring systems could be considered
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3:
        return False, "Not enough complex ring systems detected"

    return True, "Contains key features of lipopolysaccharides: trisaccharide units with long fatty acid chains"