"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:22990 pyrroline
Any organic heteromonocyclic compound with a structure based on a dihydropyrrole.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is an organic heteromonocyclic compound with a structure based on a dihydropyrrole.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a single ring with 5 atoms, 4 carbon and 1 nitrogen
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No ring found"
    rings = ring_info.AtomRings()
    pyrroline_ring = None
    for ring in rings:
        if len(ring) == 5:
            ring_atoms = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
            if ring_atoms.count(6) == 4 and ring_atoms.count(7) == 1:
                pyrroline_ring = ring
                break
    if pyrroline_ring is None:
        return False, "No pyrroline ring found"

    # Check for double bond in the ring
    bond_types = [mol.GetBondBetweenAtoms(pyrroline_ring[i], pyrroline_ring[(i+1)%5]).GetBondType() for i in range(5)]
    if Chem.rdchem.BondType.DOUBLE not in bond_types:
        return False, "No double bond found in the ring"

    # Check for exocyclic substituents
    allowed_substituents = [Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP2]
    for atom_idx in pyrroline_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in pyrroline_ring:
                if neighbor.GetHybridization() not in allowed_substituents:
                    return False, "Substituent is not allowed for pyrrolines"

    # Check for disallowed functional groups or heteroatoms
    disallowed_patterns = [Chem.MolFromSmarts("[O;!R]"), Chem.MolFromSmarts("[S;!R]")]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Molecule contains disallowed functional group or heteroatom"

    # Handle tautomeric forms and resonance structures
    tautomers = Chem.MolFromSmiles(smiles, Chem.SANITIZE_PROPERTIES, True)
    if tautomers is None:
        return False, "Invalid SMILES string"
    for tautomer in tautomers:
        if is_pyrroline(Chem.MolToSmiles(tautomer))[0]:
            return True, "Molecule is a valid pyrroline tautomer or resonance structure"

    return True, "Molecule contains a 5-membered ring with 4 carbon and 1 nitrogen atoms, a double bond, and allowed substituents"