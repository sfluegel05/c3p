"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""

from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate group (-SO3-) attached to an alkane chain.
    The sulfur is connected to a sp3-hybridized carbon atom (aliphatic, not in a ring),
    and the sulfur has three oxygen atoms: two double-bonded oxygens and one single-bonded
    oxygen with a negative charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found = False
    # Iterate over all sulfur atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Check that sulfur has degree 4 (connected to 4 atoms)
            if atom.GetDegree() != 4:
                continue
            # Sulfur should not be in a ring
            if atom.IsInRing():
                continue
            neighbors = atom.GetNeighbors()
            oxygen_count = 0
            negative_oxygen = False
            carbon_count = 0
            for neighbor in neighbors:
                # Check for oxygen neighbors
                if neighbor.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    # Check bond type
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        oxygen_count += 1
                    elif bond.GetBondType() == Chem.BondType.SINGLE:
                        # Check for negatively charged oxygen
                        if neighbor.GetFormalCharge() == -1:
                            negative_oxygen = True
                            oxygen_count += 1
                # Check for carbon neighbor
                elif neighbor.GetAtomicNum() == 6:
                    # Carbon should be sp3-hybridized and not in a ring
                    if neighbor.GetHybridization() != Chem.HybridizationType.SP3:
                        continue
                    if neighbor.IsInRing():
                        continue
                    carbon_count += 1
            # Check if sulfur is connected to exactly one carbon and three oxygens
            if carbon_count == 1 and oxygen_count == 3 and negative_oxygen:
                found = True
                break  # Found alkanesulfonate oxoanion group

    if found:
        return True, "Contains alkanesulfonate oxoanion group"
    else:
        return False, "Does not contain alkanesulfonate oxoanion group"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'alkanesulfonate oxoanion',
        'definition': 'An alkanesulfonate in which the carbon at position 1 is attached to R, which can represent hydrogens, a carbon chain, or other groups.',
        'parents': []
    }
}