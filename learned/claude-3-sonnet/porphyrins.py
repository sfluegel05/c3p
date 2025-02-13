"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:35626 porphyrins

Porphyrins are natural pigments containing a fundamental skeleton of four pyrrole nuclei
united through the alpha-positions by four methine groups to form a macrocyclic structure.
"""

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors

def is_porphyrin(smiles):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for tetrapyrrole macrocycle
    rings = mol.GetRingInfo().AtomRings()
    porphyrin_ring = None
    for ring in rings:
        if len(ring) == 16:  # Typical porphyrin macrocycle size
            porphyrin_ring = ring
            break
    if not porphyrin_ring:
        return False, "No tetrapyrrole macrocycle found"

    # Check if macrocycle has 4 pyrrole rings
    pyrrole_count = 0
    for atom_idx in porphyrin_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic() and atom.GetAtomicNum() == 7:  # Aromatic nitrogen
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:  # Connected to 2 carbons
                carbon_neighbors = [mol.GetAtomWithIdx(n.GetIdx()).GetAtomicNum() == 6 for n in neighbors]
                if all(carbon_neighbors):
                    pyrrole_count += 1
    if pyrrole_count != 4:
        return False, f"Found {pyrrole_count} pyrrole rings in macrocycle, expected 4"

    # Check for 4 methine bridges
    methine_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            atom1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
            atom2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
            if atom1.GetIsAromatic() and atom2.GetIsAromatic():
                if atom1.GetIdx() in porphyrin_ring and atom2.GetIdx() in porphyrin_ring:
                    methine_count += 1
    if methine_count != 4:
        return False, f"Found {methine_count} methine bridges, expected 4"

    # Check for conjugated system
    conjugated = True
    for atom_idx in porphyrin_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.GetIsAromatic():
            conjugated = False
            break
    if not conjugated:
        return False, "Macrocycle is not fully conjugated"

    # Check for planarity
    if not rdMolDescriptors.CalcPMI(mol, force=True) < 0.1:
        return False, "Macrocycle is not planar"

    # Additional checks for common porphyrin features
    has_substituents = any(atom.GetDegree() > 3 for atom in mol.GetAtoms() if atom.GetIdx() in porphyrin_ring)
    has_metal_center = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())

    reason = "Contains a tetrapyrrole macrocycle with 4 pyrrole rings connected by 4 methine bridges, forming a conjugated planar system"
    if has_substituents:
        reason += ", with substituents present"
    if has_metal_center:
        reason += ", and a metal center"

    return True, reason