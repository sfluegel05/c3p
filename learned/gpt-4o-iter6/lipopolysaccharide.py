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

    # Define a more complex potential polysaccharide pattern (broader sugar units)
    trisaccharide_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](O)CO)[C@H](O[C@H]1O[C@H](O)[C@@H]([C@H]1O)CO)C(O)C[C@H](O)C=O")
    if not mol.HasSubstructMatch(trisaccharide_pattern):
        return False, "No complex trisaccharide-like unit found"

    # Extended fatty acid pattern - enhanced to accommodate varying chain lengths and modifications
    lipid_chain_pattern = Chem.MolFromSmarts("C(=O)O[C@H]([CH2]10)[CH2]([CH2]1)[CH2]([CH2]1)[CH2]([CH2]1)[CH2]([CH2]1)[CH2]([CH2]1)C(=O)O")
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No adequate long fatty acid chain detected"

    # Increased complexity check for polysaccharide richness
    side_chain_pattern = Chem.MolFromSmarts("O[C@H]([C@H](O)C(O)C)[C@H](O)CO")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No sufficient oligosaccharide side chains detected"
    
    # Evaluate overall molecular complexity
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 5:
        return False, "Insufficient complex ring systems indicative of lipopolysaccharide structure"

    # Probe for presence of sugars/aromatic systems to hint polysaccharide involvement
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable_bonds < 8:
        return False, "Insufficient complexity in rotatable bonds"

    # Ensure the presence of numerous oxygen atoms indicating polysaccharide presence
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 12:
        return False, "Lack of oxygen atoms, not typical of a polysaccharide-rich lipopolysaccharide"

    return True, "Contains key structural features of lipopolysaccharides: complex polysaccharide units, side chains, and fatty acid chains"