"""
Classifies: CHEBI:26167 polar amino acid
"""
#!/usr/bin/env python
"""
Classifies: Polar Amino Acid
Definition: Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
This implementation uses a SMARTS pattern to identify a canonical amino acid backbone,
i.e. a nitrogen attached to an α‐carbon that bears a side chain and a carboxyl group.
Then, it recursively traverses the side chain (starting from the designated branch atom)
to check for polar groups (oxygen, nitrogen, or –SH).
If more than one amino acid backbone is found in the molecule we assume it is a peptide
or a larger molecule and do not classify it as a single amino acid.
"""

from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines whether a given molecule (as a SMILES string) corresponds to a polar amino acid.
    It first identifies the amino acid backbone using a SMARTS pattern with labelled atoms.
    Then, it extracts the side chain (i.e. the branch on the α‐carbon not involved in the backbone)
    to check for polar functional groups (oxygen, nitrogen, or –SH).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a polar amino acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS to capture the amino acid backbone.
    # Labeling: [N:1][C:2]([*:3])C(=O)[O:4]
    #   - [N:1]: the amino group;
    #   - [C:2]: the α‐carbon;
    #   - [*:3]: the substituent (side chain) attached to the α‐carbon;
    #   - C(=O)[O:4]: the carboxyl group.
    aa_smarts = "[N:1][C:2]([*:3])C(=O)[O:4]"
    patt = Chem.MolFromSmarts(aa_smarts)
    
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return False, "Amino acid backbone (N-αC-C(=O)O) not detected"
    
    # If there is more than one match, it likely indicates multiple amino acid units (i.e. a peptide)
    if len(matches) > 1:
        return False, "Multiple amino acid backbones detected; molecule appears to be a peptide or larger compound"
    
    # Only one match was found. Unpack the matching atom indices.
    # The SMARTS labeling gives: idx1 = nitrogen, idx2 = α-carbon, idx_sidechain = branch starting atom, idx4 = carboxyl oxygen
    idx_n, idx_alpha, idx_side, idx_o = matches[0]
    
    # For clarity, if the supposed side chain atom is hydrogen (or not heavy), we consider that there is no functional side chain.
    side_atom = mol.GetAtomWithIdx(idx_side)
    if side_atom.GetAtomicNum() < 6:
        # In the rare case the side chain is just a hydrogen (glycine), glycine is often considered non-polar.
        return False, "Side chain is hydrogen (glycine), thus not polar by our definition"
    
    # Extract the side chain subgraph by doing a depth-first search starting at the side chain atom,
    # but do not traverse back into the backbone (i.e. exclude the α-carbon).
    sidechain_atoms = set()
    stack = [idx_side]
    while stack:
        cur_idx = stack.pop()
        if cur_idx in sidechain_atoms:
            continue
        sidechain_atoms.add(cur_idx)
        cur_atom = mol.GetAtomWithIdx(cur_idx)
        for nb in cur_atom.GetNeighbors():
            # do not traverse back into the backbone (α-carbon)
            if nb.GetIdx() == idx_alpha:
                continue
            if nb.GetIdx() not in sidechain_atoms:
                stack.append(nb.GetIdx())
                
    # Analyze the side chain for polar functional groups.
    polar_found = False
    polar_features = []
    for aidx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(aidx)
        anum = atom.GetAtomicNum()
        # Oxygen or Nitrogen are intrinsically polar.
        if anum == 7 or anum == 8:
            polar_found = True
            polar_features.append(atom.GetSymbol())
        elif anum == 16:
            # For sulfur, require that at least one hydrogen is attached (–SH group)
            if atom.GetTotalNumHs() > 0:
                polar_found = True
                polar_features.append("SH")
    
    if polar_found:
        return True, f"Side chain contains polar feature(s): {', '.join(polar_features)}"
    else:
        return False, "Side chain does not contain a polar functional group capable of hydrogen bonding"

# Example usage (can be removed or kept for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "NC(CCC(N)=O)C(O)=O",       # glutamine, polar
        "NC(CO)C(O)=O",            # serine, polar
        "N[C@H](CS)C(O)=O",         # D-cysteine, polar due to –SH
        "N[C@@H](Cc1c[nH]cn1)C(O)=O",# L-histidine, polar (N in side chain)
        "NCCCC[C@@H](N)C(O)=O",     # D-lysine, polar (side chain N)
        "OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]",  # L-arginine-d7 (should be detected)
        "O(C(=O)[C@@](N([2H])[2H])(C(C(O[2H])=O)([2H])[2H])[2H])[2H]"         # L-aspartic acid-d7 (should be detected)
    ]
    for s in test_smiles:
        result, reason = is_polar_amino_acid(s)
        print(f"SMILES: {s}\nIs polar amino acid? {result} ({reason})\n")