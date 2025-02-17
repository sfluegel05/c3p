"""
Classifies: CHEBI:26167 polar amino acid
"""
#!/usr/bin/env python
"""
Classifies: Polar Amino Acid
Definition: Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
This program detects a free (non‐peptide) amino acid backbone by searching for an α‐carbon that
has exactly three heavy atom neighbors: one nitrogen (the amino group), one carboxyl carbon (which has
at least one double‐bonded oxygen) and one side‐chain branch. Glycine (side chain = H) is not considered polar,
and if multiple amino acid backbones are found the molecule is assumed to be a peptide or larger compound.
Then, the side‐chain is DFS–traversed to search for polar atoms.
"""

from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines whether the input molecule (given by a SMILES string) is a polar amino acid.
    A polar amino acid is defined as a free amino acid (not part of a peptide) 
    whose side chain has at least one polar atom (N or O, or S with at least one H).
    
    Criteria:
       1. Identify an α–carbon that has exactly three heavy-atom neighbors.
          (Heavy atoms: non‐hydrogen; note: glycine, having only two, is excluded.)
       2. The three neighbors must include:
            a. An amino group: an atom with atomic number 7.
            b. A carboxyl group: a carbon (atomic number 6) that in turn is bonded (via a double bond)
               to at least one oxygen (atomic number 8).
            c. A side-chain branch: any other heavy atom.
       3. If more than one such α–carbon is detected, return False as the molecule is likely a peptide.
       4. Traverse the side-chain (avoiding the α–carbon) and search for polar atoms:
            - Oxygen (8) or nitrogen (7).
            - Sulfur (16) but only if it has at least one attached hydrogen.
    
    Args:
         smiles (str): SMILES string representation of the molecule.
         
    Returns:
         (bool, str): Tuple of a boolean (True if it is a polar amino acid, False otherwise)
                      and a string with the reason for classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    candidates = []
    # Iterate over atoms to find candidate α–carbons (C, with three heavy-atom neighbors)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # only carbons can be α–carbon
            continue
        # Get heavy-atom neighbors (hydrogens are implicit)
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # For a (non‑glycine) amino acid, the α–carbon should have exactly three heavy neighbors.
        if len(neighbors) != 3:
            continue

        amino_found = None
        carboxyl_found = None
        sidechain_found = None
        # Examine the neighbors
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 7:
                # Found amino nitrogen candidate
                amino_found = nbr
            elif nbr.GetAtomicNum() == 6:
                # Could be the carboxyl carbon candidate
                # Check that this carbon is bonded to at least one oxygen by a double bond.
                oxy_double = False
                for nbr2 in nbr.GetNeighbors():
                    # Do not inspect the bond going back to our candidate α–carbon.
                    if nbr2.GetIdx() == atom.GetIdx():
                        continue
                    if nbr2.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            oxy_double = True
                            break
                if oxy_double:
                    carboxyl_found = nbr
                else:
                    # If not carboxyl, then treat it as side chain
                    sidechain_found = nbr
            else:
                # Any other heavy atom candidate is taken as side chain.
                sidechain_found = nbr

        # Must have one amino, one carboxyl group; if side chain is missing then it is glycine.
        if amino_found is None or carboxyl_found is None or sidechain_found is None:
            continue
        # More than one candidate α–carbon found (should be exactly one for a free amino acid)
        candidates.append({
            'alpha': atom,
            'amino': amino_found,
            'carboxyl': carboxyl_found,
            'sidechain': sidechain_found
        })
    
    if len(candidates) == 0:
        return False, "No valid free amino acid backbone (non-glycine) found"
    if len(candidates) > 1:
        return False, "Multiple amino acid backbones detected; molecule appears to be a peptide or larger compound"
    
    # Use the single candidate found
    backbone = candidates[0]
    alpha = backbone['alpha']
    side_start = backbone['sidechain']
    
    # Now traverse the side chain starting from side_start; avoid going back to the α–carbon
    side_atoms = set()
    stack = [side_start.GetIdx()]
    while stack:
        curr_idx = stack.pop()
        if curr_idx in side_atoms:
            continue
        side_atoms.add(curr_idx)
        curr_atom = mol.GetAtomWithIdx(curr_idx)
        for nbr in curr_atom.GetNeighbors():
            # Do not traverse back to the α–carbon or outside (if already visited)
            if nbr.GetIdx() == alpha.GetIdx():
                continue
            if nbr.GetIdx() not in side_atoms:
                stack.append(nbr.GetIdx())
    
    # Check for polar features in the side chain.
    polar_features = []
    for idx in side_atoms:
        atom = mol.GetAtomWithIdx(idx)
        anum = atom.GetAtomicNum()
        if anum in (7, 8):  # nitrogen or oxygen
            polar_features.append(atom.GetSymbol())
        elif anum == 16:  # sulfur; require it has at least one hydrogen (–SH)
            # GetTotalNumHs counts explicit + implicit hydrogens.
            if atom.GetTotalNumHs() > 0:
                polar_features.append("SH")
    
    if polar_features:
        return True, "Side chain contains polar feature(s): " + ", ".join(polar_features)
    else:
        return False, "Side chain does not contain polar functional group(s) capable of hydrogen bonding"

# Example test cases from the provided list.
if __name__ == "__main__":
    test_cases = [
        ("NC(CCC(N)=O)C(O)=O", "glutamine"),
        ("NC(CO)C(O)=O", "serine"),
        ("N[C@H](CS)C(O)=O", "D-cysteine"),
        ("N[C@@H](Cc1c[nH]cn1)C(O)=O", "L-histidine"),
        ("N[C@@H](CCC(N)=O)C(O)=O", "L-glutamine"),
        ("OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]", "L-arginine-d7"),
        ("NCCCC[C@@H](N)C(O)=O", "D-lysine"),
        ("N[C@H](CC(O)=O)C(O)=O", "D-aspartic acid"),
        ("NC(Cc1ccc(O)cc1)C(O)=O", "tyrosine"),
        ("N[C@H](CO)C(O)=O", "D-serine"),
        ("N[C@H](Cc1c[nH]c2ccccc12)C(O)=O", "tryptophan"),
        ("N[C@@H](CCCNC(N)=N)C(O)=O", "L-arginine")
    ]
    
    for smi, name in test_cases:
        result, reason = is_polar_amino_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result} ({reason})\n")