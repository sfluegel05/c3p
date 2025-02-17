"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected by peptide linkages.
This function has been improved: instead of simply counting any C(=O)N match and then all alpha–carbons,
we now examine each candidate peptide (amide) bond and check for the presence of an alpha–carbon candidate
on each side. For linear dipeptides a single valid peptide bond (connecting two “α–carbon” candidates) is expected.
For cyclic dipeptides (diketopiperazines), we allow a 6‐membered ring that contains two such valid peptide bonds.
"""

from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dipeptide.
    For linear dipeptides, the molecule should contain one peptide bond connecting two residues.
    For cyclic dipeptides (diketopiperazines), expect a 6‐membered ring containing two peptide bonds.
    
    To improve on our previous approach we now try to identify a “valid peptide bond”
    by checking that (a) for the carbonyl carbon (C in C(=O)N) there is at least one neighboring
    carbon (other than the peptide nitrogen) that is not terminal (heuristically an α–carbon candidate),
    and (b) for the peptide nitrogen there is at least one neighboring carbon (other than the carbonyl carbon)
    that is not terminal.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if recognized as a dipeptide, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS pattern for an amide bond (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_matches = list(mol.GetSubstructMatches(amide_pattern))
    if not peptide_matches:
        return False, "No peptide (amide) bond pattern detected."

    # Helper functions:
    def has_alpha_candidate_on_carbonyl(c_idx, n_idx):
        """
        For the carbonyl carbon (at index c_idx), check if at least one neighbor (which is carbon)
        other than the peptide nitrogen (n_idx) is present and is not terminal (degree > 1).
        This is a rough indicator for the alpha–carbon in the first amino acid residue.
        """
        atom = mol.GetAtomWithIdx(c_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != n_idx:
                # Require that the candidate carbon has degree > 1 (i.e. attached to more than a single atom).
                if nbr.GetDegree() > 1:
                    return True
        return False

    def has_alpha_candidate_on_nitrogen(n_idx, c_idx):
        """
        For the peptide nitrogen (at index n_idx), check if at least one neighbor (which is carbon)
        other than the carbonyl carbon (c_idx) is present and is not terminal.
        This is a rough indicator for the alpha–carbon in the second amino acid residue.
        """
        atom = mol.GetAtomWithIdx(n_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c_idx:
                if nbr.GetDegree() > 1:
                    return True
        return False

    # Now iterate over the amide matches and select those that show valid connectivity on both sides.
    valid_peptide_bonds = []  # each item is a tuple (c_idx, n_idx)
    for match in peptide_matches:
        # In the SMARTS "C(=O)N", match[0] is the carbonyl carbon and match[1] is the peptide nitrogen.
        c_idx, n_idx = match[0], match[1]
        if has_alpha_candidate_on_carbonyl(c_idx, n_idx) and has_alpha_candidate_on_nitrogen(n_idx, c_idx):
            valid_peptide_bonds.append((c_idx, n_idx))
    
    # Now we decide:
    # (1) Linear dipeptide: one valid peptide bond.
    if len(valid_peptide_bonds) == 1:
        return True, "Dipeptide detected: linear peptide bond connecting two α–carbon candidates."
    
    # (2) Allow cyclic dipeptides (diketopiperazines): look for a 6-membered ring that contains 2 valid peptide bonds.
    ring_info = mol.GetRingInfo()
    candidate_ring_found = False
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            # Count valid peptide bonds that are completely inside this ring.
            count_in_ring = 0
            for (c_idx, n_idx) in valid_peptide_bonds:
                if c_idx in ring and n_idx in ring:
                    count_in_ring += 1
            if count_in_ring == 2:
                candidate_ring_found = True
                break
    if candidate_ring_found:
        return True, "Cyclic dipeptide (diketopiperazine) detected: 6-membered ring with 2 peptide bonds identified."

    # If we have more than one valid peptide bond but no proper 6-membered ring detected,
    # the molecule likely has additional amide(s) (e.g. from protecting groups or additional modifications).
    if len(valid_peptide_bonds) > 1:
        return False, f"Found {len(valid_peptide_bonds)} valid peptide bond candidates; ambiguous for a dipeptide."
    
    return False, "No valid peptide bond connecting two α–carbon candidates was found."

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        # Expected dipeptides:
        "S(=O)(=N)(CC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)O",  # gamma-glutamylmethionine sulfoximine
        "SCC(NC(=O)C(N)C(O)C)C(O)=O",  # Threoninyl-Cysteine
        "CC(C)[C@H](NC(=O)[C@H](C)N)C(O)=O",  # Ala-Val
        # Some that previously caused challenges:
        "NCC(N1[C@@H](CCC1)C(=O)NC=2C=CC3=C(C2)C=CC=C3)=O",  # glycyl-L-proline 2-naphthylamide (was false negative before)
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H](SC[C@H]([NH3+])C(=O)NCC([O-])=O)[C@@H](O)CCCC([O-])=O"  # leukotriene D4, a false positive example
    ]
    for sm in test_smiles:
        result, reason = is_dipeptide(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")