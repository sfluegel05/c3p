"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: glycosaminoglycan
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Find all ring information
    ri = mol.GetRingInfo()
    ring_atom_sets = ri.AtomRings()
    sugar_ring_count = 0
    amino_sugar_count = 0

    # Define SMARTS patterns
    # General pattern for sugar ring (5 or 6-membered ring with one oxygen)
    sugar_ring_smarts = Chem.MolFromSmarts("[#6,#8]1[#6,#8][#6,#8][#6,#8][#6,#8][#6,#8]1")
    amino_group_smarts = Chem.MolFromSmarts("[#6]-[NX3;H2,H1]")  # Carbon attached to NH2 or NH

    # Set to keep track of sugar rings
    sugar_ring_atoms = []

    # Identify sugar rings
    for ring_atoms in ring_atom_sets:
        ring_size = len(ring_atoms)
        if ring_size not in (5, 6):
            continue
        # Extract the substructure of the ring
        ring_mol = Chem.PathToSubmol(mol, ring_atoms)
        # Check if ring matches sugar ring pattern
        if ring_mol.HasSubstructMatch(sugar_ring_smarts):
            sugar_ring_count += 1
            sugar_ring_atoms.extend(ring_atoms)

    if sugar_ring_count < 2:
        return False, f"Only found {sugar_ring_count} sugar ring(s), not a polysaccharide"

    # Identify amino sugars
    amino_sugar_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() in sugar_ring_atoms:
            if atom.GetAtomicNum() == 6:  # Carbon atom
                # Check for amino group attached to carbon
                if atom.HasSubstructMatch(amino_group_smarts):
                    amino_sugar_atoms.add(atom.GetIdx())

    # Count amino sugars
    amino_sugar_count = len(amino_sugar_atoms)

    if amino_sugar_count == 0:
        return False, "No amino sugars found in the sugar rings"

    proportion = amino_sugar_count / sugar_ring_count
    if proportion < 0.3:
        return False, f"Amino sugars make up {proportion:.0%} of sugar units, not substantial"

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