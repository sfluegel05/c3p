"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected via peptide linkages.
This implementation accepts both linear dipeptides (one peptide bond connecting two residues)
and cyclic dipeptides (diketopiperazines – a 6‐membered ring with 2 peptide bonds).
"""

from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if the given SMILES string corresponds to a dipeptide.
    
    For a linear dipeptide, one valid peptide bond (C(=O)-N) should connect two residues.
    For cyclic dipeptides (diketopiperazines), a 6‐membered ring containing two peptide bonds 
    is acceptable.
    
    It uses a SMARTS pattern with atom mapping to identify an amide bond.
    (Note: The SMARTS "[C:1](=O)[N:2]" yields a match tuple with three atoms where the oxygen 
    is included even though it is not mapped. Hence, we take match[0] as the carbonyl carbon
    and match[-1] as the peptide nitrogen.)
    
    For each candidate peptide bond, the function checks that the carbonyl carbon has a neighboring 
    carbon (other than the nitrogen) that can act as an α–carbon for the first residue, and 
    the peptide nitrogen has a neighboring carbon (other than the carbonyl carbon) to act as the 
    α–carbon for the second residue.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a dipeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # SMARTS for a peptide bond. Note: despite atom mapping to carbon (1) and nitrogen (2),
    # the pattern includes the oxygen, so a match tuple typically has 3 atoms.
    amide_pattern = Chem.MolFromSmarts("[C:1](=O)[N:2]")
    peptide_matches = list(mol.GetSubstructMatches(amide_pattern))
    if not peptide_matches:
        return False, "No amide (peptide) bond pattern detected."
    
    # Helper functions to check for neighboring carbon atoms (α–carbon candidates)
    def has_alpha_candidate_on_carbonyl(c_idx, n_idx):
        """
        For the carbonyl carbon at c_idx, check for a neighboring carbon (other than the peptide nitrogen at n_idx).
        """
        atom = mol.GetAtomWithIdx(c_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != n_idx:
                return True
        return False

    def has_alpha_candidate_on_nitrogen(n_idx, c_idx):
        """
        For the peptide nitrogen at n_idx, check for a neighboring carbon (other than the carbonyl carbon at c_idx).
        """
        atom = mol.GetAtomWithIdx(n_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c_idx:
                return True
        return False

    valid_peptide_bonds = []
    # Iterate over all peptide bond matches
    for match in peptide_matches:
        # Instead of unpacking into 2 values, pick the first and last atoms from the match tuple.
        if len(match) < 2:
            continue
        c_idx = match[0]       # carbonyl carbon
        n_idx = match[-1]      # peptide nitrogen
        if has_alpha_candidate_on_carbonyl(c_idx, n_idx) and has_alpha_candidate_on_nitrogen(n_idx, c_idx):
            valid_peptide_bonds.append((c_idx, n_idx))
    
    # Decision:
    # (1) For a linear dipeptide, exactly one valid peptide bond should be present.
    if len(valid_peptide_bonds) == 1:
        return True, "Dipeptide detected: a single peptide bond connects two amino acid residues."
    
    # (2) Check for cyclic dipeptide: a 6-membered ring (diketopiperazine) containing two peptide bonds.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            count_in_ring = 0
            for (c_idx, n_idx) in valid_peptide_bonds:
                if c_idx in ring and n_idx in ring:
                    count_in_ring += 1
            if count_in_ring == 2:
                return True, "Cyclic dipeptide (diketopiperazine) detected: 6-membered ring with 2 peptide bonds."
    
    # Ambiguous or unexpected cases
    if len(valid_peptide_bonds) > 1:
        return False, f"Found {len(valid_peptide_bonds)} valid peptide bond candidates; ambiguous for a dipeptide."
    
    return False, "No valid peptide bond connecting two α–carbon candidates was found."

# Example usage for testing purposes:
if __name__ == "__main__":
    test_smiles = [
        # gamma-glutamylmethionine sulfoximine (should be rejected as a dipeptide)
        "S(=O)(=N)(CC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)O",
        # glycyl-L-proline 2-naphthylamide
        "NCC(N1[C@@H](CCC1)C(=O)NC=2C=CC3=C(C2)C=CC=C3)=O",
        # Threoninyl-Cysteine
        "SCC(NC(=O)C(N)C(O)C)C(O)=O",
        # Ala-Val
        "CC(C)[C@H](NC(=O)[C@H](C)N)C(O)=O",
        # cyclic dipeptide example: cyclo(L-Leu-L-Pro)
        "C1C[C@@]2([H])C(N[C@H](C(N2C1)=O)CC(C)C)=O"
    ]
    for sm in test_smiles:
        result, reason = is_dipeptide(sm)
        print(f"SMILES: {sm}\nResult: {result}\nReason: {reason}\n")