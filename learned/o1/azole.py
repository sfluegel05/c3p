"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    Azoles can also be part of fused ring systems and may contain other heteroatoms like oxygen or sulfur.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an azole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    ring_atoms = ri.AtomRings()
    found_azole = False

    # Iterate over each ring
    for ring in ring_atoms:
        if len(ring) != 5:
            continue  # Skip rings that are not five-membered

        # Check if all atoms in the ring are aromatic
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_aromatic:
            continue  # Skip non-aromatic rings

        # Check if the ring contains at least one nitrogen atom
        contains_nitrogen = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'N' for idx in ring)
        if not contains_nitrogen:
            continue  # Skip rings without nitrogen

        # Check if all atoms in ring are C, N, O, S
        acceptable_atoms = {'C', 'N', 'O', 'S'}
        ring_atom_symbols = {mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring}
        if not ring_atom_symbols.issubset(acceptable_atoms):
            continue  # Ring contains atoms other than C, N, O, S

        # Found an azole ring (including fused rings)
        found_azole = True
        break

    if found_azole:
        return True, "Contains a five-membered aromatic ring with nitrogen (azole)"
    else:
        return False, "No five-membered aromatic ring with nitrogen found (not an azole)"