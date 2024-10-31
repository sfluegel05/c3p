from rdkit import Chem
from rdkit.Chem import AllChem

def is_naphthyridine(smiles: str):
    """
    Determines if a molecule is a naphthyridine (bicyclic compound with naphthalene skeleton 
    where two carbons are replaced by nitrogens).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a naphthyridine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 10 atoms in the rings
    rings = mol.GetRingInfo()
    ring_atoms = set()
    for ring in rings.AtomRings():
        ring_atoms.update(ring)
    if len(ring_atoms) != 10:
        return False, "Must contain exactly 10 atoms in fused ring system"

    # Check for two fused 6-membered rings
    if len(rings.AtomRings()) != 2 or not all(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "Must contain exactly two fused 6-membered rings"

    # Check aromaticity
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms):
        return False, "All ring atoms must be aromatic"

    # Count nitrogens in ring system
    n_count = 0
    n_positions = []
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'N':
            n_count += 1
            n_positions.append(atom_idx)

    if n_count != 2:
        return False, "Must contain exactly two nitrogen atoms in ring system"

    # Check remaining atoms are carbons
    if not all(mol.GetAtomWithIdx(i).GetSymbol() in ['C', 'N'] for i in ring_atoms):
        return False, "Ring system must contain only carbon and nitrogen atoms"

    # Check total number of atoms matches ring atoms
    if mol.GetNumAtoms() != len(ring_atoms):
        return False, "Molecule must consist only of the ring system"

    # All checks passed - this is a valid naphthyridine
    return True, "Valid naphthyridine structure"
# Pr=1.0
# Recall=1.0