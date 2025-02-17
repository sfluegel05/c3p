"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
A beta-lactam is defined as a lactam in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
This implementation first selects the largest fragment, then iterates over all 4-membered rings found in
the molecule. For each ring, it checks if the ring contains a nitrogen and a carbon that is part of a carbonyl
(i.e. has a double bond to an oxygen). It then ensures that this nitrogen and carbon are directly connected.
If so, the molecule is classified as a beta-lactam.
"""

from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    This function processes the SMILES, selects the largest fragment (to avoid misclassification with salts),
    and then examines each 4-membered ring to identify if it contains an amide bond in the ring (i.e. a bond
    between a nitrogen and a carbon that is double-bonded to oxygen).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-lactam, False otherwise.
        str: A reason supporting the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove salts or multiple fragments by selecting the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments could be parsed from SMILES"
    main_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Get the ring information from the main fragment.
    ring_info = main_frag.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # A tuple of tuples containing atom indices for each ring.
    
    # Loop over all rings
    for ring in atom_rings:
        if len(ring) != 4:
            continue  # Only interested in 4-membered rings.
        
        # Collect atoms in the ring.
        ring_atoms = [main_frag.GetAtomWithIdx(idx) for idx in ring]
        
        # Identify nitrogen atoms in the ring.
        nitrogens = [atom for atom in ring_atoms if atom.GetAtomicNum() == 7]
        if not nitrogens:
            continue  # A beta-lactam must contain a nitrogen.
        
        # Identify any carbon in the ring that is carbonyl.
        # We define a carbonyl carbon as one that is a carbon (atomic number 6) with at least one double bond to oxygen.
        carbonyl_carbons = []
        for atom in ring_atoms:
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbors for a double-bonded oxygen.
            for nbr in atom.GetNeighbors():
                # Skip if neighbor not oxygen.
                if nbr.GetAtomicNum() != 8:
                    continue
                bond = main_frag.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond and bond.GetBondTypeAsDouble() == 2.0:  # BondType DOUBLE: we check if it is a double bond.
                    carbonyl_carbons.append(atom)
                    break  # Found a carbonyl for this atom.
        
        if not carbonyl_carbons:
            continue  # No carbonyl group in this ring
        
        # Now check if any carbonyl carbon is directly bonded to any nitrogen in the ring.
        for c in carbonyl_carbons:
            for n in nitrogens:
                bond = main_frag.GetBondBetweenAtoms(c.GetIdx(), n.GetIdx())
                if bond:
                    # We have found a nitrogen-carbon bond inside a 4-membered ring where the carbon is carbonyl.
                    return True, ("Beta-lactam ring detected: found 4-membered ring with a carbonyl carbon " 
                                  "and an adjacent nitrogen representing an amide bond")
    
    # If we get here, no 4-membered ring meets the beta-lactam criteria.
    return False, "No beta-lactam ring found (4-membered ring with the required amide bond was not detected)"

# Example usage:
if __name__ == "__main__":
    # Test on several examples:
    test_smiles = [
        "O=C1CCN1",  # azetidin-2-one, simplest beta-lactam
        "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O",  # 7beta-aminodeacetoxycephalosporanic acid
        "C[C@H]1[C@@H](C(=O)N1S(=O)(=O)O)NC(=O)C(=NOC(C)(C)C(=O)O)C2=CSC(=N2)N"  # A complex beta-lactam example
    ]
    
    for smi in test_smiles:
        result, reason = is_beta_lactam(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("------")