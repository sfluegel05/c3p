"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
Definition: Any molecule that contains two amino‐acid residues connected via peptide linkages.
This function accepts both linear dipeptides (one peptide bond connecting two residues)
and cyclic dipeptides (diketopiperazines – a 6‐membered ring with 2 peptide bonds).
"""

from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if the given SMILES string corresponds to a dipeptide.
    For a linear dipeptide, one valid peptide bond (–C(=O)–N–) should connect two residues.
    For cyclic dipeptides (diketopiperazines), a 6‐membered ring containing two such peptide bonds is acceptable.

    This implementation uses an explicit SMARTS pattern with atom mapping to identify amide bonds.
    Then for each candidate amide bond it checks, on the carbonyl and amide sides respectively,
    that there is at least one neighboring carbon (aside from the linking partner)
    that can be heuristically considered as an “α–carbon” of an amino acid residue.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a dipeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Use an explicit SMARTS pattern for an amide bond with atom mapping.
    # Here atom 1 is the carbonyl carbon and atom 2 is the linked nitrogen.
    amide_pattern = Chem.MolFromSmarts("[C:1](=O)[N:2]")
    peptide_matches = list(mol.GetSubstructMatches(amide_pattern))
    if not peptide_matches:
        return False, "No amide (peptide) bond pattern detected."
    
    # Helper functions that check for an α–carbon candidate.
    # Instead of requiring a high connectivity, we simply require that
    # beside the amide bond partner, there is at least one carbon neighbor.
    def has_alpha_candidate_on_carbonyl(c_idx, n_idx):
        """
        For the carbonyl carbon (c_idx), verify it has a neighboring carbon
        (other than the peptide nitrogen at n_idx) that could be the α–carbon of the first residue.
        """
        atom = mol.GetAtomWithIdx(c_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != n_idx:
                return True
        return False

    def has_alpha_candidate_on_nitrogen(n_idx, c_idx):
        """
        For the peptide nitrogen (n_idx), verify it has a neighboring carbon
        (other than the carbonyl carbon at c_idx) that could be the α–carbon of the second residue.
        """
        atom = mol.GetAtomWithIdx(n_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c_idx:
                return True
        return False

    valid_peptide_bonds = []  # Will store tuples (carbonyl_idx, nitrogen_idx) that pass our tests.
    for match in peptide_matches:
        c_idx, n_idx = match  # match ordering given by SMARTS: atom mapping 1 then 2.
        if has_alpha_candidate_on_carbonyl(c_idx, n_idx) and has_alpha_candidate_on_nitrogen(n_idx, c_idx):
            valid_peptide_bonds.append((c_idx, n_idx))
    
    # Decision stage:
    # (1) A linear dipeptide is expected to have exactly one valid peptide bond.
    if len(valid_peptide_bonds) == 1:
        return True, "Dipeptide detected: linear peptide bond connecting two amino acid residues."
    
    # (2) Cyclic dipeptides (diketopiperazines) typically have a 6‐membered ring containing two peptide bonds.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            count_in_ring = 0
            for (c_idx, n_idx) in valid_peptide_bonds:
                if c_idx in ring and n_idx in ring:
                    count_in_ring += 1
            if count_in_ring == 2:
                return True, "Cyclic dipeptide (diketopiperazine) detected: 6-membered ring with 2 peptide bonds."
    
    # If there is more than one valid peptide bond but no 6‐membered ring was found, the molecule may be ambiguous.
    if len(valid_peptide_bonds) > 1:
        return False, f"Found {len(valid_peptide_bonds)} valid peptide bond candidates; ambiguous for a dipeptide."
    
    return False, "No valid peptide bond connecting two α–carbon candidates was found."

# Example usage for testing purposes:
if __name__ == "__main__":
    # List of some test SMILES strings for dipeptides (both linear and cyclic examples).
    test_smiles = [
        # gamma-glutamylmethionine sulfoximine
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