"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
Definition: An N-sulfonylurea is defined as a urea (–NH–C(=O)–NH–) in which one of the hydrogen atoms on a urea nitrogen
is replaced by a sulfonyl group (–S(=O)(=O)–). This motif is common in many herbicides and in some type 2 diabetes drugs.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    
    Our strategy is as follows:
    
    1. We search for a clear urea motif: a fragment with the pattern "N-C(=O)-N".
    2. For each urea fragment found, we check both urea nitrogens (the ones attached to the carbonyl)
       to see if one is directly substituted with a sulfonyl group.
    3. To call a substituent “sulfonyl”, the atom attached must be sulfur and that sulfur, in turn,
       must have at least two double-bonded oxygens (typical for a sulfonyl functional group).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an N-sulfonylurea, False otherwise.
        str: Reason for the classification.
    """
    
    # Try to get RDKit molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a basic urea SMARTS pattern: "N-C(=O)-N"
    urea_smarts = "N-C(=O)-N"
    urea_query = Chem.MolFromSmarts(urea_smarts)
    if urea_query is None:
        return None, None  # Should not occur
    
    # Find urea substructure matches.
    urea_matches = mol.GetSubstructMatches(urea_query)
    
    # Helper function: Given an atom (a candidate urea nitrogen), check if it is bound
    # to a sulfur that qualifies as a sulfonyl group.
    def has_sulfonyl_substituent(n_atom: rdchem.Atom):
        # Consider all neighbors except the carbonyl carbon (we expect one
        # connection from the urea motif already).
        for nbr in n_atom.GetNeighbors():
            # Skip if neighbor is the carbonyl carbon (with atomic number 6)
            if nbr.GetAtomicNum() == 6:
                # (Could try to check if that C is in a C(=O) group but for our pattern the urea C is the only C attached)
                continue
            # Look for a sulfur atom.
            if nbr.GetAtomicNum() == 16:
                # Check that this sulfur looks like part of a sulfonyl group:
                # It should have at least two bonds to oxygen with bond order DOUBLE.
                n_double_oxygen = 0
                for s_nbr in nbr.GetNeighbors():
                    # Only count oxygen (atomic number 8) that is double-bonded.
                    if s_nbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), s_nbr.GetIdx())
                        if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            n_double_oxygen += 1
                if n_double_oxygen >= 2:
                    # Found a valid sulfonyl substituent.
                    return True
        return False

    # Iterate over each urea match.
    # urea_matches returns tuples: (index_n1, index_C, index_n2)
    for match in urea_matches:
        left_n = mol.GetAtomWithIdx(match[0])
        right_n = mol.GetAtomWithIdx(match[2])
        
        # Check left nitrogen: if it has a substituent S (other than the carbonyl carbon) that qualifies as sulfonyl.
        if has_sulfonyl_substituent(left_n):
            return True, "Molecule contains a valid N-sulfonylurea motif: urea fragment with left N substituted by sulfonyl group."
        # Check right nitrogen:
        if has_sulfonyl_substituent(right_n):
            return True, "Molecule contains a valid N-sulfonylurea motif: urea fragment with right N substituted by sulfonyl group."
    
    # If no urea fragment had a valid sulfonyl substitution.
    return False, "N-sulfonylurea moiety not found in the molecule."

# Example usage (remove when using as a module):
if __name__ == "__main__":
    # For instance, test using glyburide SMILES:
    test_smiles = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"
    result, reason = is_N_sulfonylurea(test_smiles)
    print("Result:", result)
    print("Reason:", reason)