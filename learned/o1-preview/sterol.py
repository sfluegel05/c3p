"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: sterol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least 4 rings
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 4:
        return False, f"Contains only {len(ssr)} rings, requires at least 4"

    # Identify fused ring systems
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    fused_rings = []

    # Build adjacency list of rings
    ring_adjacency = {i: set() for i in range(len(atom_rings))}
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            if set(atom_rings[i]) & set(atom_rings[j]):
                ring_adjacency[i].add(j)
                ring_adjacency[j].add(i)

    # Find fused ring systems
    visited = set()
    for i in range(len(atom_rings)):
        if i in visited:
            continue
        stack = [i]
        fused_system = set()
        while stack:
            ring_idx = stack.pop()
            if ring_idx not in visited:
                visited.add(ring_idx)
                fused_system.add(ring_idx)
                stack.extend(ring_adjacency[ring_idx] - visited)
        if len(fused_system) >=4:
            fused_rings.append(fused_system)

    if not fused_rings:
        return False, "Does not contain fused ring system with at least 4 rings"

    # Check for hydroxy group attached to ring atom
    has_hydroxy = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1:
            # Oxygen atom with single bond (likely hydroxy group)
            neighbor = atom.GetNeighbors()[0]
            if neighbor.IsInRing():
                has_hydroxy = True
                break
    if not has_hydroxy:
        return False, "Does not have hydroxy group attached to ring"

    # Optional: Check if molecular framework is similar to cholestan-3-ol
    # Using Murcko Scaffold
    scaffold = Chem.MolToSmiles(rdMolDescriptors.Get_scaffold_mol(mol))
    cholestan_scaffold = Chem.MolToSmiles(rdMolDescriptors.Get_scaffold_mol(Chem.MolFromSmiles("C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@H]3CC=C4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C")))
    if scaffold != cholestan_scaffold:
        pass  # Not strictly necessary, sterols may have variations

    return True, "Contains fused ring system with at least 4 rings and hydroxy group attached to ring"

__metadata__ = {
    'chemical_class': {
        'name': 'sterol',
        'definition': 'Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol (additional carbon atoms may be present in the side chain).'
    }
}