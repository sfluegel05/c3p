"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: any organic heteromonocyclic compound based on a dihydropyrrole, i.e. pyrroline.
A pyrroline is defined here as a 5-membered ring containing exactly one nitrogen atom and exactly one double bond.
"""

from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule contains a pyrroline core (dihydropyrrole) based on its SMILES string.
    The pyrroline core is recognized as a 5-membered ring that has exactly one nitrogen and exactly one double bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a pyrroline moiety is found, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the set of rings from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Loop over each ring in the molecule
    for ring in ring_info:
        # Consider only 5-membered rings
        if len(ring) != 5:
            continue

        # Count the number of nitrogen atoms in the ring
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue  # Not the desired heterocycle

        # Now, count the number of double bonds among bonds wholly within this ring.
        double_bond_count = 0
        # Iterate over all bonds in the molecule
        for bond in mol.GetBonds():
            # Check if both atoms of the bond are in the ring
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring and a2 in ring:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_bond_count += 1

        # For dihydropyrrole (pyrroline) we expect exactly one double bond in the ring.
        if double_bond_count == 1:
            return True, "Found a 5-membered ring with one nitrogen and one double bond (dihydropyrrole core)."
    
    return False, "No 5-membered dihydropyrrole ring with one nitrogen and one double bond found."