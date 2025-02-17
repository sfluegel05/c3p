"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone – alpha,beta-unsaturated ketones
Definition: An enone is an alpha,beta-unsaturated ketone of general formula 
            R(1)R(2)C=CR(3)-C(=O)R(4) with R(4) ≠ H. In other words, 
            the motif must contain a C=C bond directly conjugated to a ketone 
            (and not an aldehyde).
            
This implementation first locates every carbonyl (C=O) group in the molecule.
It then checks that (a) the carbon bearing the C=O is a ketone (i.e. has exactly 
two heavy substituents, excluding the double-bonded oxygen) and (b) one of these 
substituents (the “beta” carbon) is attached by a single bond to the carbonyl carbon 
and further is involved in a C=C bond (the “alpha” to “beta” double bond).
If a valid enone fragment is found, the function returns True with a reason that includes 
the atom indices for the alpha carbon, beta carbon, and the carbonyl carbon.
Otherwise, it returns False with an explanation.
"""

from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (an alpha,beta-unsaturated ketone)
    based on its SMILES string.
    
    This function searches for the enone motif in two steps:
      1. Identify carbonyl carbons (C=O) in the molecule.
      2. For each carbonyl carbon, verify that it is a ketone 
         (i.e. it has exactly two heavy neighbors aside from the double-bonded oxygen)
         and that one of these neighbors (the beta carbon) is connected via a SINGLE bond
         and itself participates in a double bond (C=C) to an alpha carbon.
         This provides the required C=C–C(=O) fragment.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a valid enone motif is found, otherwise False.
        str: A message explaining the classification result.
    """
    
    # Parse molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms to look for a carbonyl group (C=O)
    for atom in mol.GetAtoms():
        # Check that the atom is carbon and has a double bond to oxygen
        if atom.GetSymbol() != "C":
            continue
        # Find oxygen neighbors that are double-bonded (the carbonyl oxygen)
        oxy_neighbors = []
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetSymbol() == "O" and bond.GetBondType() == Chem.BondType.DOUBLE:
                oxy_neighbors.append(nbr)
        if not oxy_neighbors:
            continue  # not a carbonyl
        
        # Ensure that this carbonyl carbon is a ketone and not an aldehyde.
        # For a ketone the carbon should be attached to exactly two heavy atoms,
        # excluding the oxygen of the C=O.
        heavy_nbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetSymbol() != "O"]
        if len(heavy_nbrs) != 2:
            # Either an aldehyde (only one heavy neighbor) or over‐substituted.
            continue
        
        # Now, check if one of the heavy neighbors (beta carbon) is connected
        # via a SINGLE bond and participates in a C=C double bond.
        for beta in heavy_nbrs:
            bond_to_beta = mol.GetBondBetweenAtoms(atom.GetIdx(), beta.GetIdx())
            # We require the bond between beta and carbonyl is a SINGLE bond.
            if bond_to_beta.GetBondType() != Chem.BondType.SINGLE:
                continue
            # For the beta atom, look for a neighbor (other than the carbonyl)
            # that forms a double bond (C=C) with beta.
            for nbr in beta.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue  # Skip the carbonyl carbon.
                bond_beta = mol.GetBondBetweenAtoms(beta.GetIdx(), nbr.GetIdx())
                if nbr.GetSymbol() == "C" and bond_beta.GetBondType() == Chem.BondType.DOUBLE:
                    # We have found a motif: alpha (nbr) = beta (beta) - carbonyl (atom)
                    # As an extra check, ensure that the carbonyl’s other heavy substituent (R(4))
                    # is not a hydrogen (i.e. it really is a ketone).
                    # Since we already required exactly two heavy neighbors, the other substituent is:
                    other = [n for n in heavy_nbrs if n.GetIdx() != beta.GetIdx()][0]
                    # (We assume that if it is a ketone, this substituent is not hydrogen.)
                    
                    reason = (
                        "Contains a conjugated enone motif (C=C–C(=O)) found at atom indices: "
                        "alpha (idx {}), beta (idx {}), and carbonyl carbon (idx {}).".format(
                            nbr.GetIdx(), beta.GetIdx(), atom.GetIdx()
                        )
                    )
                    return True, reason
    
    return False, ("Does not contain a conjugated enone motif "
                   "(no valid C=C–C(=O) fragment with a ketone carbonyl was found)")

# For testing purposes, you can uncomment the lines below and test several SMILES strings:
# test_smiles = [
#     "C=CC(=O)C",                # methyl vinyl ketone derivative, a minimal enone
#     "CC1=CC(=O)CC(C)(C)C1",      # isophorone, a cyclic enone
#     "O=C(C)C=C",                # an enal (aldehyde) – should not be classified as enone
#     "COc1cc2c(cc1C)C(=CC(=O)[C@]2(C)O)C(C)C",  # Lacinilene C 7-methyl ether (FP in earlier run)
# ]
#
# for sm in test_smiles:
#     result, msg = is_enone(sm)
#     print("SMILES:", sm)
#     print("Result:", result)
#     print("Reason:", msg, "\n")