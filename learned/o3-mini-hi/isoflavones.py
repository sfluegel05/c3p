"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: isoflavones – isoflavonoids defined as any isoflavonoid with a 
3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and substituted derivatives.

Our strategy:
  1. Use a SMARTS pattern to detect a benzopyran-4-one (chromen-4-one) fused ring system.
  2. Then, for atoms in that core, search for a substituent that forms part of an independent 
     six-membered aromatic ring (a phenyl group) that is not fused to the core.
     
If both criteria are met the molecule is classified as an isoflavone.
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone (3-aryl-1-benzopyran-4-one)
    based on its SMILES string.
    
    The algorithm:
      1. Parse the molecule.
      2. Look for a benzopyran-4-one (chromen-4-one) core using the SMARTS:
            c1ccc2c(c1)oc(=O)c(c2)
         This pattern tolerates additional substituents on the rings.
      3. For each match (set of core atoms), check each atom for adjacent atoms not in 
         the core. If one such neighbor is a carbon, is aromatic, and belongs to a clean 
         six‐membered aromatic ring (i.e. a phenyl ring not overlapping with the core),
         we accept the molecule as an isoflavone.
         
    Args:
         smiles (str): SMILES representation of the molecule.
         
    Returns:
         bool: True if the molecule is classified as an isoflavone.
         str: An explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS to match a benzopyran-4-one fused system.
    # This pattern is intended to capture the chromen-4-one core.
    core_smarts = "c1ccc2c(c1)oc(=O)c(c2)"
    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        return False, "Failed to compile SMARTS for benzopyran-4-one core"

    core_matches = mol.GetSubstructMatches(core)
    if not core_matches:
        return False, "Molecule does not contain a benzopyran-4-one core"
    
    # Get ring information for the whole molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of rings; each ring is a tuple of atom indices
    
    # Now, for each match of the benzopyran core, search for an external phenyl group.
    for match in core_matches:
        core_atom_set = set(match)
        # For each atom in the core, examine substituents (neighbors not in the core)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors that are not in the core fragment.
                if nbr.GetIdx() in core_atom_set:
                    continue
                # Needs to be a carbon and aromatic (a candidate for a phenyl substituent)
                if nbr.GetAtomicNum() != 6 or not nbr.GetIsAromatic():
                    continue
                # Check each ring that the neighbor is a member of.
                for ring in rings:
                    # The neighbor must be in a six-membered ring.
                    if nbr.GetIdx() in ring and len(ring) == 6:
                        # Ensure that the ring is entirely aromatic.
                        if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
                            # And ensure that none of the ring atoms are part of the core.
                            if core_atom_set.isdisjoint(ring):
                                return True, "Molecule contains a benzopyran-4-one core with an attached phenyl ring (3-aryl substitution)."
    return False, "Molecule does not contain a suitable 3-aryl substitution on the benzopyran-4-one core."

# Example testing:
if __name__ == "__main__":
    # A few test examples including one known isoflavone (formononetin) and one molecule that should not pass.
    test_examples = {
        "formononetin": "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",
        "Non-isoflavone example": "COc1c(C)c(O)c2oc(=O)c3ccc(O)cc3c12",
    }
    for name, smi in test_examples.items():
        result, reason = is_isoflavones(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}, Reason: {reason}\n")