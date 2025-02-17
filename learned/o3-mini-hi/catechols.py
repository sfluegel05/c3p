"""
Classifies: CHEBI:33566 catechols
"""
#!/usr/bin/env python3
"""
Classifies: Catechols (Any compound containing an o-diphenol component)

A catechol is defined as a molecule containing a benzene ring (6-membered aromatic ring)
with two adjacent hydroxyl (-OH) substituents. This implementation scans all aromatic rings,
and for each six-membered ring it looks for two adjacent atoms that are substituted with –OH.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol (contains an o-diphenol component) based on its SMILES string.
    A catechol is identified as any benzene ring (six-membered aromatic ring) that bears two adjacent -OH substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a catechol moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Make all hydrogens explicit to better detect -OH groups.
    mol = Chem.AddHs(mol)

    # Obtain information about all rings in the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # returns tuples of atom indices
    
    # Function to check if an atom in a ring has a hydroxyl (-OH) group attached.
    def has_hydroxyl(atom, ring_indices):
        # Look at neighbors that are not part of the ring (i.e. substituents)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_indices:
                # Check if neighbor is oxygen.
                if neighbor.GetAtomicNum() == 8:
                    # Check if this oxygen has at least one hydrogen bonded.
                    # With explicit hydrogens added, hydrogens appear as atoms with atomic number 1.
                    for sub in neighbor.GetNeighbors():
                        if sub.GetAtomicNum() == 1:
                            return True
        return False

    # Loop over all rings
    for ring in atom_rings:
        # We are interested only in six-membered rings.
        if len(ring) != 6:
            continue

        # Check that every atom in the ring is aromatic (benzene ring)
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # For each atom in the ring, record whether it carries an -OH substituent
        hydroxyl_flags = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            hydroxyl_flags.append(has_hydroxyl(atom, ring))
        
        # Since the ring is cyclic, check each pair of adjacent atoms.
        n = len(hydroxyl_flags)
        for i in range(n):
            if hydroxyl_flags[i] and hydroxyl_flags[(i+1) % n]:
                return True, "Contains o-diphenol (catechol) moiety on an aromatic ring"
    
    return False, "No adjacent hydroxyl groups on an aromatic ring found"

# Example testing when running standalone:
# examples = [
#     ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C", "oleuropein"),
#     ("C1=CC(=C(C(=C1O)O)[N+]([O-])=O)C", "4-methyl-3-nitrocatechol"),
#     ("C1(=CC=C(C=C1O)[N+]([O-])=O)O", "4-nitrocatechol"),
#     ("C=1(C=CC(=C(C1)O)O)/C=C/C(OCC)=O", "ethyl trans-caffeate"),
#     ("S(OC1=C(O)C=C([C@@H](O)CN)C=C1)(O)(=O)=O", "norepinephrine sulfate"),
# ]
#
# for smi, name in examples:
#     result, reason = is_catechols(smi)
#     print(f"NAME: {name} -> {result}: {reason}")