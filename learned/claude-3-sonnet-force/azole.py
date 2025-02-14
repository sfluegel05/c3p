"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:34545 azole
An azole is any monocyclic or polycyclic heterocycle consisting of a five-membered ring containing nitrogen.
Azoles can also contain one or more other non-carbon atoms, such as nitrogen, sulfur or oxygen.
"""

from rdkit import Chem
from rdkit.Chem import rdchem, rdMolDescriptors

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
    rings = mol.GetRingInfo().AtomRings()
    azole_rings = [ring for ring in rings if len(ring) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring)]

    if not azole_rings:
        return False, "No 5-membered ring containing nitrogen found"

    # Check if any of the azole rings are part of a fused ring system (polycyclic)
    polycyclic_azole_rings = []
    for ring in azole_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        fused_rings = rdchem.FindFusedRings(mol, ring_atoms)
        if len(fused_rings) > 1:
            polycyclic_azole_rings.append(ring)

    # Check for allowed heteroatoms (N, O, S)
    allowed_heteroatoms = [7, 8, 16]
    for ring in polycyclic_azole_rings or azole_rings:
        ring_heteroatoms = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() not in [1, 6]]
        if any(heteroatom not in allowed_heteroatoms for heteroatom in ring_heteroatoms):
            return False, "Ring contains disallowed heteroatoms"

    # Additional checks to ensure the azole ring is part of a valid azole moiety
    azole_pattern = Chem.MolFromSmarts("[*:1]1[*:2][*:3][*:4][*:5]1")
    azole_matches = mol.GetSubstructMatches(azole_pattern)
    for match in azole_matches:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        if any(atom.GetAtomicNum() == 7 for atom in ring_atoms):
            if any(atom.GetTotalDegree() > 3 for atom in ring_atoms):
                continue  # Avoid quaternary nitrogen or other unusual structures
            if all(atom.GetTotalDegree() == 2 for atom in ring_atoms):
                continue  # Avoid aromatic stabilized rings
            if sum(atom.GetImplicitValence() for atom in ring_atoms) > 8:
                continue  # Avoid highly unsaturated or strained rings
            if rdMolDescriptors.CalcNumAromaticRings(mol) == 1:
                continue  # Avoid monocyclic aromatic rings (e.g., pyrrole)
            return True, "Molecule contains a monocyclic or polycyclic azole moiety"

    return False, "No valid azole moiety found"