"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as a compound with a fully conjugated cyclic dione structure
    derived from aromatic compounds by conversion of an even number of -CH= groups
    into -C(=O)- groups with any necessary rearrangement of double bonds (polycyclic and heterocyclic analogues are included).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    if not ri.IsInitialized():
        return False, "No ring structures found in molecule"

    atom_rings = ri.AtomRings()
    if not atom_rings:
        return False, "No ring structures found in molecule"

    # Iterate over each ring in the molecule
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        conjugated = True
        ketone_carbons = []

        # Check if the ring is fully conjugated
        for i, atom in enumerate(ring_atoms):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%len(ring)])
            if not bond.GetIsConjugated():
                conjugated = False
                break
        if not conjugated:
            continue

        # Count the number of ketone groups (C=O) attached to carbons in the ring
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 6:  # Carbon atom
                for nbr in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone_carbons.append(atom.GetIdx())
                        break

        # Check if there are at least two ketone groups (even number)
        if len(ketone_carbons) >= 2 and len(ketone_carbons) % 2 == 0:
            return True, "Ring with fully conjugated cyclic dione structure found"

    return False, "Does not contain a quinone core structure"

__metadata__ = {
    'chemical_class': {
        'name': 'quinone',
        'definition': 'Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones, derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups with any necessary rearrangement of double bonds (polycyclic and heterocyclic analogues are included).',
    },
    'examples': [
        'O=C1C=CC(=O)C=CC1=O',  # Benzoquinone
        'O=C1C=CC2=CC=CC=C2C1=O',  # Naphthoquinone
        'O=C1C=CC2=C1C=CC(=O)C=C2',  # Anthraquinone
    ],
    'version': '2.1',
    'author': 'Assistant',
}