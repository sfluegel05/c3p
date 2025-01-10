"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain a complex mixture of specific trisaccharide units,
    oligosaccharide structures, and long fatty acid chains.

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

    # Pattern for complex trisaccharide unit (two heptose units and octulosonic acid)
    trisaccharide_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O[C@H]2O[C@H](O)[C@@H]([C@H]2O)CO)O[C@@H]([C@@H]1O)CO)[C@H](O)[C@@H](O)C(=O)CO") # Example pattern; refine as needed
    if not mol.HasSubstructMatch(trisaccharide_pattern):
        return False, "No specific trisaccharide unit found"

    # Pattern for fatty acid chains (here using 3-hydroxytetradecanoic acid as an example)
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC(=O)O") # Example pattern; adjust for hydroxyl groups
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No long fatty acid chains detected"

    # Pattern for polysaccharide side chains
    oligosaccharide_pattern = Chem.MolFromSmarts("C(CO)O[C@H]([C@H](O)CO)O") # Example pattern; extend to cover various oligosaccharides
    if not mol.HasSubstructMatch(oligosaccharide_pattern):
        return False, "No sufficient oligosaccharide side chains detected"

    # Count the number of rings as an indicator of structural complexity
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 5:
        return False, "Insufficient complex ring systems indicative of lipopolysaccharide structure"

    # Ensure the molecule has enough rotatable bonds suggesting polysaccharide structuring
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable_bonds < 10:
        return False, "Insufficient complexity in rotatable bonds"

    # Check for an abundance of oxygen atoms indicating polysaccharide presence
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 15:
        return False, "Insufficient oxygen atoms, not characteristic of a polysaccharide-rich lipopolysaccharide"

    return True, "Contains key structural features of lipopolysaccharides: complex polysaccharide units, oligosaccharide side chains, and fatty acid chains"