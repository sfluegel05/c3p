"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
Definition: An N-sulfonylurea is defined as a urea moiety (–NH–C(=O)–NH–)
in which one of the two urea nitrogens is directly substituted by a sulfonyl group (–S(=O)(=O)–).
This motif is common to several herbicides and type 2 diabetes drugs.
Our approach is:
    1. Examine each carbonyl carbon that is double-bonded to an O.
    2. Require that carbon to be part of a urea-like fragment (connected to exactly two N atoms).
    3. For each of the two N atoms (which should be neutral), check if it is substituted
       (apart from the urea C) by a sulfur.
    4. Then verify that the candidate sulfur is “sulfonyl” – that is, it has at least two double bonds to oxygen.
If such a fragment is found, we call the molecule an N-sulfonylurea.
Note: This method was refined after analysis of mis‐classified (false positive or false negative)
examples; extra filters (e.g. checking the formal charge) help avoid some cases (like anionic nitrogen).
"""

from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    
    Our strategy is to (1) search for a urea-like fragment: a carbonyl carbon (C=O)
       bearing exactly two nitrogen substituents, and then (2) for each of those N atoms,
       see whether (apart from the urea carbon) it is directly substituted by a sulfur 
       atom that carries at least two oxygen atoms by double bond.
       
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a valid N-sulfonylurea fragment is found, False otherwise.
        str: A reason (or message) describing the result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Iterate over all atoms looking for a carbonyl carbon (C=O)
    for atom in mol.GetAtoms():
        # Look at carbon atoms only
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check that the carbon has at least one double-bonded oxygen (C=O)
        oxy_double = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    oxy_double.append(nbr)
        if not oxy_double:
            continue  # not a carbonyl
        
        # Now check that the carbon is in a urea-like fragment, i.e.
        # it is connected to exactly two nitrogen atoms.
        n_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 7]
        if len(n_neighbors) != 2:
            continue  # not the central carbon of a urea
        
        # For each urea nitrogen, check for a sulfonyl substituent
        for n in n_neighbors:
            # Skip if the nitrogen is charged (we expect a neutral urea N)
            if n.GetFormalCharge() != 0:
                continue
            # Look over neighbors of this N excluding the urea carbon
            for sub in n.GetNeighbors():
                if sub.GetIdx() == atom.GetIdx():
                    continue
                # Check if neighbor is a sulfur candidate.
                if sub.GetAtomicNum() == 16:
                    # Verify that the S is attached by a SINGLE bond to our nitrogen.
                    bond_NS = mol.GetBondBetweenAtoms(n.GetIdx(), sub.GetIdx())
                    if bond_NS is None or bond_NS.GetBondType() != Chem.BondType.SINGLE:
                        continue
                    # Count how many of the S neighbors are oxygen with a double bond.
                    dO_count = 0
                    for s_nbr in sub.GetNeighbors():
                        if s_nbr.GetAtomicNum() == 8:
                            bond_SO = mol.GetBondBetweenAtoms(sub.GetIdx(), s_nbr.GetIdx())
                            if bond_SO is not None and bond_SO.GetBondType() == Chem.BondType.DOUBLE:
                                dO_count += 1
                    # Accept S if it has at least 2 double bonds to O.
                    if dO_count >= 2:
                        return True, (
                            f"Molecule contains a valid N-sulfonylurea motif: "
                            f"urea fragment (C(=O) between N atoms) with N (atom idx {n.GetIdx()}) substituted by a sulfonyl group."
                        )
                        
    return False, "N-sulfonylurea moiety not found in the molecule."


# For testing purposes (remove or comment out this part when using as a module)
if __name__ == "__main__":
    # Try one of the true positive examples (glyburide)
    test_smiles = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"
    result, reason = is_N_sulfonylurea(test_smiles)
    print("Result:", result)
    print("Reason:", reason)