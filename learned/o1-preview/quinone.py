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
    atom_rings = ri.AtomRings()
    if not atom_rings:
        return False, "No ring structures found in molecule"

    # Define a SMARTS pattern for ring carbons double-bonded to oxygen (carbonyls)
    ketone_smarts = Chem.MolFromSmarts('[cR]=O')
    if ketone_smarts is None:
        return False, "Invalid SMARTS pattern"

    # Find all carbonyl carbons in the molecule
    ketone_matches = mol.GetSubstructMatches(ketone_smarts)
    ketone_carbon_idxs = set(match[0] for match in ketone_matches)  # Get carbon atom indices

    # Iterate over each ring in the molecule
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Check if the ring is fully conjugated
        conjugated = True
        for i in range(len(ring)):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%len(ring)])
            if not bond.GetIsConjugated():
                conjugated = False
                break
        if not conjugated:
            continue  # Skip rings that are not fully conjugated

        # Count the number of carbonyl groups in the ring
        ketone_count = sum(1 for idx in ring if idx in ketone_carbon_idxs)

        # Check if there are at least two ketone groups (even number)
        if ketone_count >= 2 and ketone_count % 2 == 0:
            return True, f"Ring with fully conjugated cyclic dione structure found with {ketone_count} ketone groups"

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
    'version': '2.2',
    'author': 'Assistant',
}