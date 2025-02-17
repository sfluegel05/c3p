"""
Classifies: CHEBI:26167 polar amino acid
"""
#!/usr/bin/env python
"""
Classifies: Polar Amino Acid
Definition: Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
This implementation identifies a free (non‐peptide) amino acid backbone via a SMARTS pattern that
accounts for neutral and ionized forms, and then extracts the side chain to search for polar functional groups.
Note: Glycine, which has only a hydrogen as the side chain, is not considered polar.
"""

from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines whether the input molecule (SMILES string) is a polar amino acid.
    
    Procedure:
      1. Use a modified SMARTS to detect the free amino acid backbone that allows for both neutral 
         (NH2) and protonated (NH3+) amino groups, as well as carboxyl in either protonated (COOH)
         or deprotonated (COO-) forms.
         The pattern used is:
             "[N:1;H2,H3+][C:2;H]([*:3])C(=O)[O:4;H0,H1,-]"
         where:
             - [N:1;H2,H3+] matches NH2 or NH3+,
             - [C:2;H] is the α‐carbon,
             - [*:3] is the branch point for the side chain, and 
             - C(=O)[O:4;H0,H1,-] matches the carboxyl group.
             
      2. Gather all matches and then extract the unique α‐carbon indices (mapped as “2”).
         If more than one unique α‐carbon is found, then the molecule likely is a peptide or larger.
         
      3. If exactly one backbone is identified, check if the side-chain branch (atom labeled “3”) is more than a hydrogen.
      
      4. Traverse the side chain (avoiding the backbone α‐carbon) and search for polar atoms:
             - Oxygen (atomic number 8) or Nitrogen (7) are polar.
             - Sulfur (16) is counted as polar only if it carries at least one hydrogen (–SH).
             
      5. Return True with an appropriate explanation if a polar feature is found; otherwise, return False.

    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        (bool, str): A boolean indicating whether the molecule is a polar amino acid and a string
                     explaining the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern that captures the free amino acid backbone.
    # It allows for both NH2 and NH3+ (the [N:1]) and carboxyl groups in different states.
    aa_smarts = "[N:1;H2,H3+][C:2;H]([*:3])C(=O)[O:4;H0,H1,-]"
    patt = Chem.MolFromSmarts(aa_smarts)
    if patt is None:
        return False, "Invalid SMARTS pattern"
    
    # Get all substructure matches for the amino acid backbone
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return False, "Amino acid backbone (free form) not detected"
    
    # Extract unique α‐carbon indices (mapped as atom "2") from the matches.
    alpha_indices = set()
    mapping_matches = []
    for match in matches:
        # match order: [N, α‐carbon, side-chain attachment, carboxyl oxygen]
        alpha_idx = match[1]
        alpha_indices.add(alpha_idx)
        mapping_matches.append(match)
    
    # If more than one unique α‐carbon is found, assume a peptide or larger compound.
    if len(alpha_indices) > 1:
        return False, "Multiple amino acid backbones detected; molecule appears to be a peptide or larger compound"

    # In the accepted case, use the first match (all share the same α‐carbon)
    match = mapping_matches[0]
    idx_n = match[0]      # amino group N
    idx_alpha = match[1]  # α‐carbon
    idx_side = match[2]   # side chain branch point
    # idx_o = match[3]   # carboxyl oxygen, not further used
    
    # Check if the side chain is just hydrogen (glycine) by checking the atomic number.
    side_atom = mol.GetAtomWithIdx(idx_side)
    if side_atom.GetAtomicNum() < 6:
        return False, "Side chain is hydrogen (glycine), thus not polar by our definition"
    
    # Traverse the side chain starting from the branch atom:
    # We perform a depth-first search and avoid going back to the α‐carbon.
    sidechain_atoms = set()
    visited = set()
    stack = [idx_side]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        sidechain_atoms.add(current)
        atom = mol.GetAtomWithIdx(current)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == idx_alpha:
                continue  # do not backtrack to the backbone
            if nbr.GetIdx() not in visited:
                stack.append(nbr.GetIdx())
                
    # Check the gathered side chain atoms for polar features.
    polar_found = False
    polar_features = []
    for aidx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(aidx)
        atomic_num = atom.GetAtomicNum()
        if atomic_num in (7, 8):  # Nitrogen or Oxygen
            polar_found = True
            polar_features.append(atom.GetSymbol())
        elif atomic_num == 16:  # Sulfur: check for -SH condition (at least one hydrogen)
            if atom.GetTotalNumHs() > 0:
                polar_found = True
                polar_features.append("SH")
    
    if polar_found:
        return True, "Side chain contains polar feature(s): " + ", ".join(polar_features)
    else:
        return False, "Side chain does not contain a polar functional group capable of hydrogen bonding"

# Example test cases (a subset from the provided list)
if __name__ == "__main__":
    test_cases = [
        ("NC(CCC(N)=O)C(O)=O", "glutamine"),
        ("NC(CO)C(O)=O", "serine"),
        ("N[C@H](CS)C(O)=O", "D-cysteine"),
        ("N[C@@H](Cc1c[nH]cn1)C(O)=O", "L-histidine"),
        ("N[C@@H](CCC(N)=O)C(O)=O", "L-glutamine"),
        ("NCCCC[C@@H](N)C(O)=O", "D-lysine"),
        ("OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]", "L-arginine-d7"),
        ("N[C@H](CC(O)=O)C(O)=O", "D-aspartic acid"),
        ("NC(Cc1ccc(O)cc1)C(O)=O", "tyrosine"),
        ("N[C@H](CO)C(O)=O", "D-serine")
    ]
    
    for smi, name in test_cases:
        result, reason = is_polar_amino_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result} ({reason})\n")