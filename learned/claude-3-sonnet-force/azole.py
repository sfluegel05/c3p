"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:34545 azole
An azole is any monocyclic or polycyclic heterocycle consisting of a five-membered ring containing nitrogen.
Azoles can also contain one or more other non-carbon atoms, such as nitrogen, sulfur or oxygen.
"""

from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a 5-membered ring containing nitrogen
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    azole_rings = []
    for ring in rings:
        if len(ring) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            azole_rings.append(ring)

    if not azole_rings:
        return False, "No 5-membered ring containing nitrogen found"

    # Check for allowed heteroatoms (N, O, S)
    allowed_heteroatoms = [7, 8, 16]
    for ring in azole_rings:
        ring_heteroatoms = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() not in [1, 6]]
        if any(heteroatom not in allowed_heteroatoms for heteroatom in ring_heteroatoms):
            return False, "Ring contains disallowed heteroatoms"

    return True, "Molecule contains a monocyclic or polycyclic heterocycle with a 5-membered ring containing nitrogen and possibly other allowed heteroatoms (O, S)"