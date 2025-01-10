"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain a complex mixture of specific trisaccharide units, oligosaccharide structures, and long fatty acid chains.

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

    # Define a possible improved trisaccharide pattern (simplified)
    # Focuses on common sugar rings and their linkages
    trisaccharide_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](O)CO)[C@@H](O)C(O)C[C@H](O)C=O")
    if not mol.HasSubstructMatch(trisaccharide_pattern):
        return False, "No complex trisaccharide-like unit found"

    # Define a fatty acid mimic pattern targeting long chains (including hydroxyl groups and length robustness)
    lipid_chain_pattern = Chem.MolFromSmarts("C(=O)OCC([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]([CH2])[CH2]C(=O)O")
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No adequate long fatty acid chain detected"

    # Attempt to capture general polysaccharide richness
    # Pattern could include or target common linkages with high connectivity
    side_chain_pattern = Chem.MolFromSmarts("O[C@H]([C@H](O)C(O)C)[C@H](O)CO")  # General polysaccharide motif
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No sufficient oligosaccharide side chains detected"
    
    # Overall heuristic rate: Monitor patterns related to complexity
    n_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:
        return False, "Insufficient complex ring systems indicative of lipopolysaccharide structure"

    return True, "Contains key structural features of lipopolysaccharides: specific trisaccharide units, complex side chains, and long fatty acid chains"