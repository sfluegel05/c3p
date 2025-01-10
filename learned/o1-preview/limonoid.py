"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a triterpenoid that is highly oxygenated and has a prototypical structure
    containing or derived from a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Check for fused ring system of three 6-membered rings and one 5-membered ring (steroid nucleus)
    ri = mol.GetRingInfo()
    ring_atoms = ri.AtomRings()
    if len(ring_atoms) < 4:
        return False, "Less than 4 rings found in molecule"

    # Count ring sizes
    ring_sizes = [len(ring) for ring in ring_atoms]
    if ring_sizes.count(6) < 3 or ring_sizes.count(5) < 1:
        return False, "Does not contain fused 6-6-6-5 ring system"

    # Step 2: Check that rings are fused properly
    # Create a map of which rings share atoms
    fused_counts = [0]*len(ring_atoms)
    for i, ring_i in enumerate(ring_atoms):
        for j, ring_j in enumerate(ring_atoms):
            if i >= j:
                continue
            if set(ring_i) & set(ring_j):
                fused_counts[i] +=1
                fused_counts[j] +=1
    # A proper steroid nucleus will have rings fused such that each ring is fused to at least two others
    if not all(count >= 2 for count in fused_counts):
        return False, "Rings are not properly fused like steroid nucleus"

    # Step 3: Check for 4,4,8-trimethyl substitutions on the ring system
    # Identify methyl groups attached to the ring system at tertiary carbons
    methyl_substitutions = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 1:  # Methyl group carbon
            neighbor = atom.GetNeighbors()[0]
            if neighbor.IsInRing():
                # Check if neighbor carbon is tertiary (attached to three carbons)
                if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 4:
                    methyl_substitutions += 1
    if methyl_substitutions < 3:
        return False, f"Only {methyl_substitutions} methyl groups attached to ring system, need at least 3"

    # Step 4: Check for furan ring attached to the ring system
    furan_pattern = Chem.MolFromSmarts('c1ccoc1')
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found in molecule"

    # Check if furan ring is connected to the ring system
    furan_atoms = set()
    for match in furan_matches:
        furan_atoms.update(match)
    ring_system_atoms = set()
    for ring in ring_atoms:
        ring_system_atoms.update(ring)
    connected = False
    for furan_atom in furan_atoms:
        atom = mol.GetAtomWithIdx(furan_atom)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in ring_system_atoms:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Furan ring is not connected to the ring system"

    # Step 5: Check if molecule is highly oxygenated (e.g., at least 6 oxygen atoms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Not highly oxygenated, contains only {o_count} oxygen atoms"

    # Step 6: Look for common functional groups in limonoids (e.g., lactones, acetates)
    lactone_pattern = Chem.MolFromSmarts('O=C1OC[*]1')  # Simple lactone pattern
    acetate_pattern = Chem.MolFromSmarts('C(=O)OC')     # Acetate ester pattern
    has_lactone = mol.HasSubstructMatch(lactone_pattern)
    has_acetate = mol.HasSubstructMatch(acetate_pattern)
    if not (has_lactone or has_acetate):
        return False, "No common limonoid functional groups (lactones or acetates) found"

    return True, "Molecule is a limonoid (contains fused 6-6-6-5 ring system with methyl substitutions and attached furan ring)"