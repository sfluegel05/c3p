"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose in our definition is a hexose (total 6 carbons and 6 oxygens)
that contains a sugar ring (5- or 6-membered with exactly one oxygen) and 
has exactly one exocyclic –CH2OH substituent attached to a ring-carbon.
That ring carbon (the “C5” center) must have CIP configuration 'R'.
"""

from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines whether a molecule is a D-hexose (with the D-configuration at C5)
    based on its SMILES string.

    Requirements:
      - Must contain exactly 6 carbons and 6 oxygens.
      - Must feature a sugar ring that is either 5- or 6-membered and contains exactly one oxygen.
      - Must contain exactly one exocyclic CH2OH substituent that is attached to a ring carbon.
        This ring carbon (candidate "C5") should have CIP configuration 'R'.

    Args:
        smiles (str): SMILES string for the molecule.

    Returns:
        (bool, str): Tuple of classification result and an explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the hydrogen counts are unambiguous.
    mol = Chem.AddHs(mol)
    # Assign stereochemistry to compute CIP codes.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Count heavy atoms: exactly 6 carbons and exactly 6 oxygens.
    atoms = list(mol.GetAtoms())
    carbons = [a for a in atoms if a.GetAtomicNum() == 6]
    oxygens = [a for a in atoms if a.GetAtomicNum() == 8]
    if len(carbons) != 6:
        return False, "Molecule does not have exactly 6 carbon atoms (not a hexose)"
    if len(oxygens) != 6:
        return False, "Molecule does not have exactly 6 oxygen atoms (modified hexose?)"
    
    # Find a candidate sugar ring.
    # We expect a 5- or 6-membered ring that contains exactly one oxygen (the ring oxygen).
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring = None
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms in the ring.
        count_ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if count_ring_oxygens == 1:
            sugar_ring = set(ring)
            break
    if sugar_ring is None:
        return False, "No appropriate sugar ring (5- or 6-membered with one oxygen) found"
    
    # Now identify candidate exocyclic CH2OH groups.
    # We are looking for a carbon atom (not in the ring) that:
    #   - Has exactly 2 (implicit+explicit) hydrogens (i.e. is a CH2 group)
    #   - Has a single oxygen attached via a single bond (i.e. the –OH)
    #   - Is attached to one of the sugar ring carbons.
    exo_candidates = []
    for atom in mol.GetAtoms():
        # Only consider carbon atoms that are NOT in the sugar ring.
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetIdx() in sugar_ring:
            continue
        # Check that this exocyclic candidate has exactly 2 hydrogens
        # (using the total hydrogen count which is set after AddHs)
        if atom.GetTotalNumHs() != 2:
            continue
        # Check that it is attached to exactly one oxygen via a single bond.
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondTypeAsDouble() == 1:
                    oxy_neighbors.append(nbr)
        if len(oxy_neighbors) != 1:
            continue
        # Check that this candidate is attached to a sugar ring carbon.
        attached_ring_atom = None
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in sugar_ring:
                attached_ring_atom = nbr
                break
        if attached_ring_atom is None:
            continue
        # Save the ring carbon that carries the exocyclic CH2OH group.
        exo_candidates.append(attached_ring_atom)
    
    if len(exo_candidates) == 0:
        return False, "No appropriate exocyclic CH2OH substituent found"
    if len(exo_candidates) > 1:
        return False, "Multiple exocyclic CH2OH substituents found (ambiguous candidate for C5)"
    
    # The single candidate ring carbon (C5) is expected to have CIP configuration 'R'.
    candidate_C5 = exo_candidates[0]
    if not candidate_C5.HasProp('_CIPCode'):
        return False, "Stereochemistry was not assigned for the candidate C5 atom"
    cip = candidate_C5.GetProp('_CIPCode')
    if cip != 'R':
        return False, f"C5 has CIP configuration {cip} (expected R for D-hexose)"
    
    return True, "Molecule is a D-hexose: contains a valid sugar ring with an exocyclic CH2OH (C5) having R configuration"

# Example usage and testing:
if __name__ == "__main__":
    test_smiles_list = [
        # Examples from the prompt (they should be accepted if they meet our criteria)
        "O1[C@H]([C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)CO",            # beta-D-idofuranose
        "O1[C@@H]([C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO",            # beta-D-altropyranose
        "[H][C@]1(O[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)CO",            # beta-D-galactofuranose
        "O1[C@@H]([C@H](O)[C@H](O)[C@H](O)[C@H]1O)CO",              # alpha-D-talopyranose
        "O1[C@@H]([C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)CO",             # beta-D-idopyranose
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",                # alpha-D-galactose
        "O1[C@@H]([C@@H](O)[C@H](O)[13C@@H](O)[13C@H]1O)CO",         # alpha-D-glucose-1,2-((13)C2)
        "[C@H]1(O)C(O)O[C@H](C=O)[C@H]([C@@H]1O)O",                 # 6-dehydro-D-glucose
        "C1([C@H]([C@@H]([C@H]([C@H](O1)CO)O)O)O)O",                # D-idopyranose
        "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O",              # beta-D-talopyranose
    ]
    for smi in test_smiles_list:
        result, reason = is_D_hexose(smi)
        print("SMILES: ", smi)
        print("Result: ", result)
        print("Reason: ", reason)
        print("-" * 50)