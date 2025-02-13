"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: ChEBI:31535 Azole
An azole is any monocyclic heteroarene consisting of a five-membered ring containing nitrogen. 
Azoles can also contain one or more other non-carbon atoms, such as nitrogen, sulfur or oxygen.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

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

    # Find all 5-membered rings containing nitrogen
    azole_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 5 and any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            azole_rings.append(set(ring))

    # Check if there is exactly one heterocyclic ring (azole)
    heterocyclic_rings = set(ring for ring in mol.GetRingInfo().AtomRings() if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring))
    if len(heterocyclic_rings) != 1 or list(heterocyclic_rings)[0] not in azole_rings:
        return False, "Molecule does not contain exactly one monocyclic azole ring"

    # Check if the azole ring contains additional heteroatoms (N, O, S)
    azole_ring = list(heterocyclic_rings)[0]
    has_additional_heteroatoms = any(mol.GetAtomWithIdx(idx).GetAtomicNum() not in [6, 7] for idx in azole_ring)

    if has_additional_heteroatoms:
        return True, "Molecule contains a monocyclic azole ring with additional heteroatoms (N, O, S)"
    else:
        return True, "Molecule contains a monocyclic azole ring without additional heteroatoms"