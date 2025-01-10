"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is a steroid with an alpha-oriented hydroxyl group at the C17 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has stereochemistry information
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Check for steroid backbone (four fused rings: three six-membered, one five-membered)
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    ring_counts = {}
    for ring in atom_rings:
        ring_size = len(ring)
        ring_counts[ring_size] = ring_counts.get(ring_size, 0) + 1

    if ring_counts.get(6, 0) < 3 or ring_counts.get(5, 0) < 1:
        return False, "Does not have the typical steroid ring structure"

    # Check for hydroxyl group at a chiral tertiary carbon (possible C17 position)
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@H](O)[C]")
    matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not matches:
        return False, "No chiral carbon with attached hydroxyl group found"

    # Check if any of these carbons could be the C17 position
    found = False
    for match in matches:
        idx_C = match[0]
        atom_C = mol.GetAtomWithIdx(idx_C)
        # Ensure the carbon is tertiary (connected to three carbons)
        neighbor_carbons = [nbr for nbr in atom_C.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(neighbor_carbons) == 3:
            # Check stereochemistry
            chirality = atom_C.GetChiralTag()
            if chirality == Chem.CHI_TETRAHEDRAL_CCW:
                found = True
                break

    if not found:
        return False, "No tertiary chiral carbon with alpha-hydroxyl group found"

    return True, "Contains steroid backbone with 17alpha-oriented hydroxyl group"