"""
Classifies: CHEBI:33566 catechols
"""
#!/usr/bin/env python3
"""
Classifies: Catechols 
Definition: Any compound containing an o-diphenol component.
A catechol is defined as an aromatic benzene ring (six-membered ring where every atom is aromatic)
bearing two adjacent substituents that are hydroxyl groups (–OH or their deprotonated/esterified form)
provided that the oxygen is directly bound to the aromatic carbon via a single bond and not further bonded to a carbon.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol moiety (an o-diphenol component)
    based on its SMILES string.
    
    The algorithm works by:
      1. Converting the SMILES to an RDKit molecule.
      2. For every bond between two aromatic carbons that are part of at least one six‐membered ring, 
         checking if both atoms have a substituent oxygen attached by a single bond.
      3. In qualifying the oxygen substituent, we require that aside from being bonded to the aromatic carbon,
         the oxygen is not further bonded to any carbon (so that –OCH3 is excluded), but may be attached to a hydrogen
         (as in –OH) or to another heteroatom (e.g. as in sulfate esters).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains an o-diphenol (catechol) moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Obtain ring information (list of tuples: each is a set of atom indices in a ring)
    rings = mol.GetRingInfo().AtomRings()

    # Helper: does a given aromatic carbon atom (in a six-membered ring) have a qualifying oxygen substituent?
    def has_qualifying_oxygen(atom):
        for neighbor in atom.GetNeighbors():
            # Exclude atoms that are part of the ring when attached via a bond in a six-membered aromatic ring.
            # We are looking for an external substituent.
            if neighbor.GetAtomicNum() != 8:
                continue
            # The bond between the aromatic atom and the oxygen must be single (exclude e.g. carbonyl oxygens)
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Now check that this oxygen is NOT further bonded to a carbon (other than the aromatic atom).
            qualifies = True
            for o_nb in neighbor.GetNeighbors():
                # skip the aromatic atom itself
                if o_nb.GetIdx() == atom.GetIdx():
                    continue
                # If the oxygen is attached to any carbon, then it is likely a methoxy (or similar), not a free hydroxyl.
                # (Exceptions such as sulfate esters are allowed since S, for example, is not carbon.)
                if o_nb.GetAtomicNum() == 6:
                    qualifies = False
                    break
            if qualifies:
                return True
        return False

    # Now iterate over all bonds in the molecule.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Consider only bonds connecting two aromatic carbons.
        if not (a1.GetIsAromatic() and a2.GetIsAromatic()):
            continue

        # Check if the two atoms a1 and a2 are part of a six‐membered ring.
        in_six_membered_ring = False
        for ring in rings:
            # Only consider rings of exactly six atoms.
            if len(ring) == 6 and (a1.GetIdx() in ring) and (a2.GetIdx() in ring):
                in_six_membered_ring = True
                break
        if not in_six_membered_ring:
            continue

        # For this aromatic bond (i.e. adjacent atoms on the ring), check for qualifying oxygen substituents.
        if has_qualifying_oxygen(a1) and has_qualifying_oxygen(a2):
            return True, "Contains o-diphenol (catechol) moiety on a six-membered aromatic ring"

    return False, "No adjacent qualifying hydroxyl substituents found on a six-membered aromatic ring"

# Example usage (for testing):
if __name__ == '__main__':
    examples = [
        ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C", "oleuropein"),
        ("C1=CC(=C(C(=C1O)O)[N+]([O-])=O)C", "4-methyl-3-nitrocatechol"),
        ("C1(=CC=C(C=C1O)[N+]([O-])=O)O", "4-nitrocatechol"),
        ("C=1(C=CC(=C(C1)O)O)/C=C/C(OCC)=O", "ethyl trans-caffeate"),
        ("S(OC1=C(O)C=C([C@@H](O)CN)C=C1)(O)(=O)=O", "norepinephrine sulfate"),
        ("O[C@H](CC\\C=C\\c1ccc(O)cc1)CCc1ccc(O)c(O)c1", "(-)-(3S)-1-(3,4-dihydroxyphenyl)-7-(4-hydroxyphenyl)-(6E)-6-hepten-3-ol"),
        ("CC(C)c1cccc(O)c1O", "3-isopropylcatechol"),
    ]
    for smi, name in examples:
        result, reason = is_catechols(smi)
        print(f"NAME: {name} -> {result}: {reason}")