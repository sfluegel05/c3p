"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
A beta-lactam is defined as a lactam (cyclic amide) in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
This implementation selects the largest fragment from the SMILES (to minimize interference from salts)
and then iterates over all 4-membered rings. For each ring, it checks:
  - the ring has exactly 4 atoms (by definition),
  - the ring has exactly one nitrogen and three carbons,
  - at least one of the carbons has a double bond to oxygen (with the oxygen not being in the ring),
  - and that carbon (the carbonyl carbon) is directly bonded to the nitrogen.
If all these conditions are met, the molecule is classified as a beta-lactam.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    The function first sanitizes the molecule, selects the largest fragment, adds explicit hydrogens,
    and then inspects each 4-membered ring. A valid beta-lactam ring must have exactly one nitrogen
    and three carbons, one of which is carbonyl (i.e. double-bonded to an oxygen that is not part of the ring)
    and is directly bonded to that nitrogen.
    
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
    
    # Remove salts/multiple fragments by selecting the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments could be parsed from the SMILES"
    
    # Choose the fragment with the most heavy atoms
    main_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Add explicit hydrogens to help with bond perceptions.
    main_frag = Chem.AddHs(main_frag)
    
    # Get ring information for the main fragment
    ring_info = main_frag.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # a tuple of tuples, each inner tuple is a set of atom indices for that ring.
    
    # Loop over all rings in the main fragment.
    for ring in atom_rings:
        if len(ring) != 4:
            continue  # We only care about 4-membered rings.
        
        # Get the atoms of this ring.
        ring_atoms = [main_frag.GetAtomWithIdx(idx) for idx in ring]
        
        # Check composition: A true beta-lactam ring should have exactly 1 nitrogen and 3 carbons.
        num_nitrogen = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        num_carbon = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        if num_nitrogen != 1 or num_carbon != 3:
            continue  # This ring does not have the correct atom composition
        
        # Identify the nitrogen atom in the ring.
        nitrogen_atom = [atom for atom in ring_atoms if atom.GetAtomicNum() == 7][0]
        
        # Look among the carbons for a carbonyl carbon.
        # A carbonyl carbon is defined here as a carbon that is double-bonded to an oxygen atom,
        # where the oxygen is not part of the 4-membered ring.
        carbonyl_carbon = None
        for atom in ring_atoms:
            if atom.GetAtomicNum() != 6:
                continue  # We are only interested in carbons.
            # Examine each bond from this carbon.
            for nbr in atom.GetNeighbors():
                # Skip if the neighbor is within the ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:  # Oxygen candidate
                    bond = main_frag.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    # Check if the bond is a double bond.
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        carbonyl_carbon = atom
                        break
            if carbonyl_carbon:
                break
        
        if not carbonyl_carbon:
            continue  # No carbonyl carbon found in this ring
        
        # Ensure the carbonyl carbon and the nitrogen are directly bonded (i.e. form the amide linkage within the ring).
        bond = main_frag.GetBondBetweenAtoms(nitrogen_atom.GetIdx(), carbonyl_carbon.GetIdx())
        if bond is not None:
            return True, ("Beta-lactam ring detected: found a 4-membered ring with one nitrogen, three carbons "
                          "(where one carbon is carbonyl via a double bond to oxygen outside the ring) and a direct "
                          "C-N bond representing the amide bond")
    
    # If no ring met our criteria, return False.
    return False, "No beta-lactam ring found (4-membered ring with the required amide bond was not detected)"

# Example usage:
if __name__ == "__main__":
    # A few example SMILES (including azetidin-2-one as the simplest beta-lactam)
    test_smiles = [
        "O=C1CCN1",   # azetidin-2-one
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O",  # 7beta-aminodeacetoxycephalosporanic acid
        "C[C@H]1[C@@H](C(=O)N1S(=O)(=O)O)NC(=O)C(=NOC(C)(C)C(=O)O)C2=CSC(=N2)N"  # complex beta-lactam derivative
    ]
    
    for smi in test_smiles:
        result, reason = is_beta_lactam(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("------")