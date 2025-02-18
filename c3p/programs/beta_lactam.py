"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
A beta-lactam is defined as a lactam (cyclic amide) in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
This implementation first selects the largest fragment, adds explicit hydrogens, and then iterates over all 4-membered rings.
For each 4-membered ring, it checks:
  - The ring has exactly four atoms.
  - The ring has exactly one nitrogen and three carbon atoms.
  - Among the carbon atoms, one has a double bond to an oxygen atom that is not part of the ring.
  - The nitrogen and that carbon (carbonyl carbon) are directly bonded by a single bond (the amide bond).
If all these conditions are met for any ring, the molecule is classified as a beta-lactam.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    The function first sanitizes the molecule, selects the largest fragment, adds explicit hydrogens,
    and then inspects each 4-membered ring. A valid beta-lactam ring must have exactly one nitrogen
    and three carbons, one of which is carbonyl (i.e. double-bonded to an oxygen outside the ring)
    and that carbon must be directly bonded by a single bond to the nitrogen.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-lactam, False otherwise.
        str: A reason supporting the classification decision.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all fragments; select the largest (to remove solvents/salts).
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments could be parsed from the SMILES"
    main_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Add explicit hydrogens to improve bond perception and handle hydrogens in amide bonds.
    main_frag = Chem.AddHs(main_frag)
    
    # Obtain ring information from the main fragment.
    ring_info = main_frag.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Search over all rings looking for a valid beta-lactam ring.
    for ring in atom_rings:
        # Only consider 4-membered rings.
        if len(ring) != 4:
            continue
        
        # Retrieve the atoms in the ring.
        ring_atoms = [main_frag.GetAtomWithIdx(idx) for idx in ring]
        
        # Count nitrogen and carbon atoms in the ring.
        num_nitrogen = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        num_carbon = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        if num_nitrogen != 1 or num_carbon != 3:
            continue
        
        # Identify the nitrogen atom of the ring.
        nitrogen_atom = next(atom for atom in ring_atoms if atom.GetAtomicNum() == 7)
        
        # Look among the carbons in the ring for a candidate carbonyl carbon.
        candidate_carbonyl = None
        for atom in ring_atoms:
            if atom.GetAtomicNum() != 6:
                continue  # only interested in carbons
            # Check neighbors of this carbon to find an oxygen making a double bond
            # Note: ensure this oxygen is not in the ring.
            oxygen_double_found = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    bond = main_frag.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        oxygen_double_found = True
                        break
            if oxygen_double_found:
                candidate_carbonyl = atom
                break
        
        if candidate_carbonyl is None:
            continue  # no carbonyl group in this ring
        
        # Check that the candidate carbonyl carbon and the nitrogen are directly bonded.
        bond = main_frag.GetBondBetweenAtoms(nitrogen_atom.GetIdx(), candidate_carbonyl.GetIdx())
        if bond is None:
            continue
        # The bond representing the amide linkage should be a single bond.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        
        # We found a ring fulfilling all criteria.
        return True, ("Beta-lactam ring detected: found a 4-membered ring with one nitrogen and three carbons, "
                      "where one carbon bears a double-bonded oxygen (outside the ring) and is directly bonded "
                      "to the nitrogen via a single bond representing the amide linkage")
    
    return False, "No beta-lactam ring found (4-membered ring with the required amide bond was not detected)"


# Example usage:
if __name__ == "__main__":
    # A small set of example SMILES including the simplest beta-lactam.
    test_smiles = [
        "O=C1CCN1",   # azetidin-2-one, simplest beta-lactam
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O",  # 7beta-aminodeacetoxycephalosporanic acid
        "C[C@H]1[C@@H](C(=O)N1S(=O)(=O)O)NC(=O)C(=NOC(C)(C)C(=O)O)C2=CSC(=N2)N"  # complex beta-lactam derivative
    ]
    
    for smi in test_smiles:
        result, reason = is_beta_lactam(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("------")