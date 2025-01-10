"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: glycosaminoglycan
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to ensure proper valence perception
    try:
        Chem.Kekulize(mol)
    except Chem.KekulizeException:
        pass

    # Detect sugar rings using a more comprehensive pattern
    # Monosaccharide rings (pyranose and furanose)
    sugar_smarts = Chem.MolFromSmarts("""
        [C;R]1([O;R][C;R][C;R][C;R][C;R]1) | 
        [C;R]1([O;R][C;R][C;R][C;R]1)
    """)
    sugar_matches = mol.GetSubstructMatches(sugar_smarts)
    sugar_rings = set()
    for match in sugar_matches:
        ring = set(match)
        sugar_rings.add(frozenset(ring))

    sugar_ring_count = len(sugar_rings)

    if sugar_ring_count < 2:
        return False, f"Only found {sugar_ring_count} sugar ring(s), not a polysaccharide"

    # Identify amino sugars
    amino_sugar_count = 0
    for ring_atoms in sugar_rings:
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                # Look for attached amino group
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                        # Check if it is an amino group (NH2 or NH)
                        num_H = neighbor.GetTotalNumHs()
                        if num_H > 0 and neighbor.GetDegree() <= 3:
                            amino_sugar_count += 1
                            break

    if amino_sugar_count == 0:
        return False, "No amino sugars found in the sugar rings"

    proportion = amino_sugar_count / sugar_ring_count
    if proportion < 0.3:
        return False, f"Amino sugars make up {proportion:.0%} of sugar units, not substantial"

    # Check for glycosidic linkages between sugar rings
    # Glycosidic bond pattern: sugar ring oxygen connected to carbon of another sugar ring
    glycosidic_bonds = 0
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetIdx() in ring_atoms and atom2.GetIdx() in ring_atoms:
            if bond.IsInRing():
                continue
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                glycosidic_bonds += 1

    if glycosidic_bonds < sugar_ring_count - 1:
        return False, "Not enough glycosidic linkages to form a polysaccharide"

    return True, f"Molecule is a polysaccharide with {amino_sugar_count} amino sugar(s) out of {sugar_ring_count} sugar rings"

__metadata__ = {
    'chemical_class': {
        'name': 'glycosaminoglycan',
        'definition': 'Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.'
    },
    'config': {
        # Configuration details can be added here
    },
    'message': None,
    'success': True,
    'error': '',
    'stdout': None
}