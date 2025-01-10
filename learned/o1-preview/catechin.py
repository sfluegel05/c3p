"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a hydroxyflavan that has a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()
    if num_rings < 3:
        return False, f"Not enough rings ({num_rings} found), catechin has 3 rings"

    # Identify the rings
    rings = ri.AtomRings()
    aromatic_rings = []
    heterocyclic_oxygen_rings = []
    for ring in rings:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        atomic_nums = set(atom.GetAtomicNum() for atom in atoms_in_ring)
        is_aromatic = all(atom.GetIsAromatic() for atom in atoms_in_ring)
        if is_aromatic:
            aromatic_rings.append(set(ring))
        elif 8 in atomic_nums:
            heterocyclic_oxygen_rings.append(set(ring))

    if len(aromatic_rings) < 2:
        return False, f"Not enough aromatic rings ({len(aromatic_rings)} found), catechin has 2 aromatic rings"
    if len(heterocyclic_oxygen_rings) < 1:
        return False, "No heterocyclic oxygen-containing ring found, catechin has an oxygen-containing heterocyclic ring"

    # Check connections between rings
    # Catechin has two aromatic rings connected via a heterocyclic ring containing oxygen
    connected = False
    for het_ring in heterocyclic_oxygen_rings:
        connections = 0
        for arom_ring in aromatic_rings:
            # Find bonds between heterocyclic ring and aromatic ring
            for atom_idx in het_ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in arom_ring:
                        connections += 1
                        break
            if connections >= 2:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Aromatic rings are not properly connected via heterocyclic ring"

    # Check for hydroxyl groups
    # At least two hydroxyl groups attached to aromatic carbons
    num_pheno_oh = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                num_pheno_oh += 1
    if num_pheno_oh < 2:
        return False, f"Not enough phenolic hydroxyl groups ({num_pheno_oh} found), catechins typically have at least two"

    # Check for hydroxyl group at position 3 (attached to saturated carbon in heterocyclic ring)
    has_3_oh = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            # Check if neighbor carbon is in heterocyclic ring and is saturated (degree 4, not aromatic)
            if neighbor.GetAtomicNum() == 6 and neighbor.IsInRing():
                if not neighbor.GetIsAromatic() and neighbor.GetDegree() == 4:
                    has_3_oh = True
                    break
    if not has_3_oh:
        return False, "Hydroxyl group at position 3 not found"

    return True, "Molecule matches catechin structure with necessary rings and hydroxyl groups"