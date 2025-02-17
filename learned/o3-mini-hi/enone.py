"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone – alpha,beta-unsaturated ketones
Definition: A molecule with an alpha,beta-unsaturated ketone of general formula
R(1)R(2)C=CR(3)-C(=O)R(4) (with R(4) not being H), that is, the carbonyl is conjugated
to a C=C double bond.
This program uses local neighborhood analysis rather than a single rigid SMARTS.
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone must have a conjugated C=C–C(=O)–C fragment:
      • A ketone carbonyl carbon (C double bonded to O) that is not an aldehyde.
      • One of the carbonyl carbon's carbon neighbors (alpha carbon) must be sp2,
        with a double bond to another carbon (the beta carbon).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as an enone, otherwise False.
        str: Reason explaining the classification.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Loop through atoms to find candidate ketone carbonyl carbons.
    for atom in mol.GetAtoms():
        # Look for carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue

        # Identify double-bonded oxygen(s)
        oxygens = []
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2:
                oxygens.append(nbr)
        if not oxygens:
            continue  # Not a carbonyl carbon if no C=O bond.

        # To be a ketone, the carbonyl carbon must be bonded to at least two carbons.
        # (An aldehyde would have one hydrogen.)
        non_oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 8]
        if len(non_oxygen_neighbors) < 2:
            continue  # Likely an aldehyde, not an enone.

        # Now check if one of these carbon neighbors is an alpha carbon
        # that is attached via a single bond and itself forms a double bond with another carbon.
        for alpha in non_oxygen_neighbors:
            bond_alpha = mol.GetBondBetweenAtoms(atom.GetIdx(), alpha.GetIdx())
            if bond_alpha.GetBondType() != Chem.BondType.SINGLE:
                continue
            if alpha.GetAtomicNum() != 6:
                continue  # Must be carbon

            # Look at neighbors of alpha (except the carbonyl carbon) for a double bond.
            for beta in alpha.GetNeighbors():
                if beta.GetIdx() == atom.GetIdx():
                    continue
                bond_alpha_beta = mol.GetBondBetweenAtoms(alpha.GetIdx(), beta.GetIdx())
                if bond_alpha_beta.GetBondType() == Chem.BondType.DOUBLE and beta.GetAtomicNum() == 6:
                    # We have: beta = carbon (double bonded to alpha),
                    # alpha = carbon (singly bonded to carbonyl),
                    # atom = carbonyl carbon (double bonded to oxygen) with at least one extra carbon neighbor.
                    return True, ("Contains a conjugated enone motif: an alpha,beta-unsaturated "
                                  "ketone via carbonyl carbon (idx {}) and alpha carbon (idx {}).".format(atom.GetIdx(), alpha.GetIdx()))
    # If no candidate was found, return False.
    return False, "Does not contain a conjugated enone motif (no C=C–C(=O)–C arrangement was found)"

# Example usage (for testing - uncomment if needed):
# test_smiles = [
#     "C=CC(=O)C",  # minimal enone (methyl vinyl ketone derivative)
#     "CC1=CC(=O)CC(C)(C)C1",  # e.g., isophorone, should match enone motif
#     "O=C(C)C=C",  # an enal (but aldehyde, should not be classified as enone)
# ]
# for sm in test_smiles:
#     result, reason = is_enone(sm)
#     print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")