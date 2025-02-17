"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
Definition: Any aldehyde in which the carbonyl group is attached to an aromatic (benzenoid) moiety.
For example, piperonal ([H]C(=O)c1ccc2OCOc2c1) and salicylaldehyde ([H]C(=O)c1ccccc1O) qualify.
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde where the carbonyl carbon (CHO group)
    is exocyclic and its sole alkyl substituent is an aromatic carbon that belongs
    to an entirely carbon six‐membered ring (a benzene or benzenoid ring).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a qualifying arenecarbaldehyde group, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Precompute benzene (or benzene-like) rings.
    # We define these as rings of size 6 with all atoms being carbon.
    ring_info = mol.GetRingInfo()
    benzene_ring_atom_sets = []
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                benzene_ring_atom_sets.append(set(ring))
    # Combine all carbon indices that are part of any benzene ring.
    benzene_atoms = set()
    for ring_set in benzene_ring_atom_sets:
        benzene_atoms |= ring_set
    
    found_aldehyde = False
    # Iterate over atoms. We look for an aldehyde carbon:
    # a carbon that is not itself aromatic, has exactly one explicit hydrogen,
    # exactly two heavy neighbors, and one of them is oxygen via a double bond.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # In a proper aldehyde the carbon should be CHO (one hydrogen)
        # and should be sp2 so it is not flagged aromatic.
        if atom.GetTotalNumHs() != 1 or atom.GetIsAromatic():
            continue
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(neighbors) != 2:
            continue
        # Identify if one neighbor is a double-bonded oxygen.
        oxygen_found = False
        aromatic_neighbor_ok = False
        other_neighbor = None
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                oxygen_found = True
            else:
                other_neighbor = nbr
        if not oxygen_found or (other_neighbor is None):
            continue
        # We have an aldehyde candidate (CHO group)
        found_aldehyde = True
        # Now check that the non-oxygen neighbor is an aromatic carbon that
        # belongs to a benzene (six-membered, all-carbon) ring.
        if other_neighbor.GetAtomicNum() == 6 and other_neighbor.GetIsAromatic():
            # Check if this neighbor belongs to one of the benzene rings.
            if other_neighbor.GetIdx() in benzene_atoms:
                return True, "Arenecarbaldehyde functional group detected (aldehyde connected to benzenoid aromatic moiety)."
    
    # If we have found at least one aldehyde group but none qualifies,
    # then report that the aldehyde is not connected to a benzene ring.
    if found_aldehyde:
        return False, "Aldehyde group present, but not connected to a benzenoid aromatic moiety."
    else:
        return False, "No aldehyde group detected."

# Example usage (uncomment to test):
# print(is_arenecarbaldehyde("[H]C(=O)c1ccc2OCOc2c1"))  # piperonal should return True
# print(is_arenecarbaldehyde("CC(=O)C"))  # aliphatic aldehyde should return False