"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone – alpha,beta-unsaturated ketones
Definition: A molecule is an enone if it contains an alpha,beta-unsaturated ketone fragment,
            which we define as a conjugated C=C–C(=O)–C unit (with the carbonyl carbon not being an aldehyde).
In this program we use local neighborhood analysis enhanced with hybridization and explicit hydrogen checks.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone (alpha,beta-unsaturated ketone) based on its SMILES string.
    
    The search is performed by identifying a candidate ketone carbonyl carbon that:
      - Is sp2 hybridized and forms a double bond to oxygen (a carbonyl)
      - Has no attached hydrogen (thus not an aldehyde)
    and then finding a neighbor (alpha carbon) that is:
      - Bound to the carbonyl via a single bond and sp2 hybridized
      - Possessing at least one neighbor (other than the carbonyl) that forms a double bond (beta carbon) 
        which is also sp2 hybridized.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if an enone motif is found, otherwise False.
        str: A reason explaining the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through atoms looking for candidate ketone carbonyl carbons.
    for atom in mol.GetAtoms():
        # Focus only on carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue

        # Check if atom is sp2 (as expected for a carbonyl group)
        if atom.GetHybridization() != rdchem.HybridizationType.SP2:
            continue

        # Look for a double-bonded oxygen to mark a carbonyl.
        oxygen_neighbors = []
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Check if bond is a double bond and neighbor is oxygen.
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                oxygen_neighbors.append(nbr)
        if not oxygen_neighbors:
            continue  # no carbonyl group here

        # To be a ketone, the carbonyl carbon should not have any explicit hydrogen (i.e. not an aldehyde).
        if atom.GetTotalNumHs() > 0:
            continue

        # Now, get the non-oxygen (heavy) neighbors.
        non_oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 8]
        if len(non_oxygen_neighbors) < 2:
            continue  # Likely an aldehyde or other non-ketone

        # Loop over neighbors to identify an alpha carbon.
        for alpha in non_oxygen_neighbors:
            # The bond from the carbonyl carbon to the alpha carbon must be a SINGLE bond.
            bond_ca = mol.GetBondBetweenAtoms(atom.GetIdx(), alpha.GetIdx())
            if bond_ca.GetBondType() != Chem.BondType.SINGLE:
                continue

            # We require the alpha carbon to also be sp2 (olefinic).
            if alpha.GetHybridization() != rdchem.HybridizationType.SP2:
                continue

            # Now search for a partner beta: a neighbor of alpha (excluding the candidate carbonyl)
            # that is attached via a DOUBLE bond and is sp2 hybridized.
            for beta in alpha.GetNeighbors():
                if beta.GetIdx() == atom.GetIdx():
                    continue
                bond_ab = mol.GetBondBetweenAtoms(alpha.GetIdx(), beta.GetIdx())
                if bond_ab.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                # Beta should be a carbon and also sp2.
                if beta.GetAtomicNum() != 6:
                    continue
                if beta.GetHybridization() != rdchem.HybridizationType.SP2:
                    continue

                # Found a candidate motif: beta = alpha,beta double bond, alpha connected by single bond to carbonyl.
                reason = ("Contains a conjugated enone motif: an alpha,beta-unsaturated ketone via carbonyl carbon (idx {}) "
                          "and alpha carbon (idx {}).".format(atom.GetIdx(), alpha.GetIdx()))
                return True, reason

    return False, "Does not contain a conjugated enone motif (no C=C–C(=O)–C arrangement was found)"

# For basic testing, you can uncomment the following lines:
# test_smiles = [
#     "C=CC(=O)C",                # minimal enone (methyl vinyl ketone derivative)
#     "CC1=CC(=O)CC(C)(C)C1",      # isophorone, should match enone motif
#     "O=C(C)C=C",                # an enal (aldehyde) – should not be classified as enone
# ]
# for sm in test_smiles:
#     result, msg = is_enone(sm)
#     print("SMILES:", sm)
#     print("Result:", result)
#     print("Reason:", msg, "\n")