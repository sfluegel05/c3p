"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
#!/usr/bin/env python
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
Definition: An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl
groups are located at positions 1 and 2 respectively.
The improved strategy is to first locate the phosphocholine headgroup (using a SMARTS fragment).
Then we anchor the glycerol backbone assuming a linear chain:
  sn-3 (CH2) attached to the phosphocholine oxygen →
  sn-2 (central carbon) attached to an acyl (ester) oxygen →
  sn-1 (terminal carbon) attached to an alkyl (ether) oxygen.
We require that the acyl branch shows a proper carbonyl (C=O) and that the alkyl branch does not.
Some flexibility (e.g. for a methyl at sn-1) is allowed.
"""

from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES.
    Args:
        smiles (str): SMILES representation of the molecule.
    Returns:
        bool: True if the molecule matches the class, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # --- STEP 1: Locate the phosphocholine headgroup ---
    # The SMARTS pattern is for a phosphocholine group. It is assumed that:
    #   [C]O P(=O)([O-])OCC[N+](C)(C)C
    # Here we expect that the oxygen that connects to the glycerol backbone is in position 1 (index 1).
    phospho_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phospho_frag = Chem.MolFromSmarts(phospho_smarts)
    phospho_matches = mol.GetSubstructMatches(phospho_frag)
    if not phospho_matches:
        return False, "Phosphocholine headgroup not found"
    
    # --- HELPER FUNCTIONS ---
    def is_acyl_branch(o_atom):
        """
        Checks if an oxygen atom is attached to a carbon that bears a carbonyl (C=O)
        indicating an ester (acyl) branch.
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                # check if there is at least one double bond from this carbon to an oxygen atom
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    def is_alkyl_branch(o_atom, min_carbons=1):
        """
        Checks if an oxygen atom leads to an alkyl (ether) branch.
        It does a DFS from the oxygen neighbor (a carbon) and counts carbons;
        if a carbonyl (C=O) is encountered in the branch, we do not consider it as an alkyl branch.
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            visited = set()
            stack = [nbr]
            c_count = 0
            while stack:
                current = stack.pop()
                if current.GetIdx() in visited:
                    continue
                visited.add(current.GetIdx())
                # if we see a carbonyl bond along this branch, then this is not a pure alkyl branch
                for b in current.GetBonds():
                    if b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = b.GetOtherAtom(current)
                        if other.GetAtomicNum() == 8:
                            return False
                # count carbon only
                if current.GetAtomicNum() == 6:
                    c_count += 1
                # if we have reached at least the minimal number, we can give a positive result
                if c_count >= min_carbons:
                    return True
                for nn in current.GetNeighbors():
                    if nn.GetAtomicNum() == 6 and nn.GetIdx() not in visited:
                        stack.append(nn)
        return False

    # --- STEP 2: Identify the glycerol backbone based on the phosphocholine linkage ---
    # Our strategy: starting from the phosphocholine match, we know that the bridging oxygen in the 
    # SMARTS (index 1) should be bound to a carbon of the glycerol (sn-3).
    for match in phospho_matches:
        # match indices (from the SMARTS) are assumed to be:
        # match[0]: carbon (preceding the bridging oxygen)
        # match[1]: oxygen bridging to phosphorus (candidate for connection to glycerol sn-3)
        # match[2]: phosphorus, etc.
        o_phospho = mol.GetAtomWithIdx(match[1])
        
        # Get the glycerol candidate attached to this oxygen.
        # It should be a carbon (sn-3) and not the phosphorus.
        sn3_candidates = [nbr for nbr in o_phospho.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != match[0]]
        if not sn3_candidates:
            continue
        for sn3 in sn3_candidates:
            # We expect sn3 to be a primary carbon (CH2 typically)
            # Next, from sn3 we try to trace to sn-2. Typically, sn-3 is connected to sn-2 (and only to o_phospho).
            sn3_neighbors = [atom for atom in sn3.GetNeighbors() if atom.GetAtomicNum() == 6 and atom.GetIdx() != o_phospho.GetIdx()]
            if len(sn3_neighbors) != 1:
                # if ambiguous, skip candidate
                continue
            sn2 = sn3_neighbors[0]
            # Now, from sn2, we expect a second backbone carbon (sn-1) besides sn3.
            sn2_neighbors = [atom for atom in sn2.GetNeighbors() if atom.GetAtomicNum() == 6 and atom.GetIdx() != sn3.GetIdx()]
            if not sn2_neighbors:
                continue
            # In many glycerol backbones there is exactly one candidate for sn1. If more, try each.
            for sn1 in sn2_neighbors:
                # --- Check substituents ---
                # (a) Check sn-2: must have an oxygen (not linking back to sn1 or sn3) that leads to an acyl branch.
                acyl_found = False
                for nbr in sn2.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [sn1.GetIdx(), sn3.GetIdx()]:
                        if is_acyl_branch(nbr):
                            acyl_found = True
                            break
                if not acyl_found:
                    continue

                # (b) Check sn-1: must have an oxygen (not linking to sn2) that leads to an alkyl branch.
                alkyl_found = False
                for nbr in sn1.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != sn2.GetIdx():
                        if is_alkyl_branch(nbr):
                            alkyl_found = True
                            break
                if not alkyl_found:
                    continue

                # (c) Additional check for sn-3: it should only be bonded to sn-2 and the phosphocholine oxygen.
                sn3_backbone = [a for a in sn3.GetNeighbors() if a.GetAtomicNum() == 6]
                if len(sn3_backbone) != 1:
                    continue

                # If we have matched the backbone in the expected order with proper branches, accept.
                return True, ("Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure with glycerol backbone "
                              "(sn-1: ether/alkyl branch, sn-2: acyl branch, sn-3: phosphocholine linkage)")
    
    return False, "No glycerol backbone with required substituents (acyl at sn-2, alkyl at sn-1, and phosphocholine at sn-3) found"


# Example usage: testing a few SMILES strings from the provided list.
if __name__ == "__main__":
    test_smiles_list = [
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC(=O)[H])COCCCCCCCCCCCCCCCC",  # 1-hexadecyl-2-formyl-sn-glycero-3-phosphocholine
        "CCCCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCC\\C=C/C\\C=C/C\\C=C/CCCCC",  # 1-hexadecyl-2-[(8Z,11Z,14Z)-eicosatrienoyl]-sn-glycero-3-phosphocholine
        "P(OC[C@@H](COCCCCCCCCCCCCCCCC)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(=O)(OCC[N+](C)(C)C)[O-]",  # 1-O-hexadecyl-2-arachidonoyl-sn-glycero-3-phosphocholine
        "C(OP(=O)(OCC[N+](C)(C)C)[O-])[C@@H](COC)OC(CCCCCCCCCCCCCCC)=O",  # 1-methyl-2-hexadecanoyl-sn-glycero-3-phosphocholine (previously false negative)
    ]
    
    for smi in test_smiles_list:
        result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")