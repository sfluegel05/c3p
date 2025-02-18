"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
#!/usr/bin/env python
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
Definition: An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl
groups are located at positions 1 and 2 respectively.
The strategy is to first locate the phosphocholine headgroup (using a SMARTS fragment),
then trace the glycerol backbone starting from the oxygen bound to phosphorus.
The glycerol backbone is expected to be a three-carbon chain where:
  • The carbon (sn-3) directly linked to the phosphocholine oxygen is CH2.
  • The neighboring carbon (sn-2) carries an oxygen that serves as an acyl (ester) branch,
    i.e. an oxygen attached to a carbon which bears a C=O.
  • The other carbon (sn-1) is bound to an oxygen that leads to a long alkyl (ether) chain.
If these conditions all hold, then the structure is accepted.
"""

from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES.
    Args:
        smiles (str): SMILES representation of the molecule
    Returns:
        bool: True if the molecule matches the class, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- STEP 1: Search for the phosphocholine headgroup ---
    # We use a SMARTS pattern that is expected to match an oxygen attached to phosphorus and to a choline fragment.
    phospho_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phospho_frag = Chem.MolFromSmarts(phospho_smarts)
    phospho_matches = mol.GetSubstructMatches(phospho_frag)
    if not phospho_matches:
        return False, "Phosphocholine headgroup not found"

    # --- HELPER FUNCTIONS ---
    def is_acyl_branch(o_atom):
        """
        Checks if an oxygen atom leads to an acyl (ester) branch.
        We require that the oxygen is attached to a carbon which is double-bonded to another oxygen.
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                # look for a double bond C=O from this carbon
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    def is_alkyl_branch(o_atom, min_carbons=5):
        """
        Checks if an oxygen atom leads to an alkyl (ether) branch.
        We perform a simple depth-first search (DFS) starting from the neighbor carbon
        to see if we can count at least 'min_carbons' carbon atoms.
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
                if current.GetAtomicNum() == 6:
                    c_count += 1
                if c_count >= min_carbons:
                    return True
                for nn in current.GetNeighbors():
                    # Avoid going back to the original oxygen atom.
                    if nn.GetIdx() not in visited and nn.GetAtomicNum() == 6:
                        stack.append(nn)
        return False

    # --- STEP 2: Identify the glycerol backbone by "walking" from phosphocholine ---
    # Our assumption: in 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, the phosphocholine attaches via an oxygen
    # to the sn-3 carbon of a glycerol backbone.
    for match in phospho_matches:
        # The SMARTS "COP(=O)([O-])OCC[N+](C)(C)C" yields a match where:
        #   match[0] : a carbon (from the "C" right before the bridging oxygen)
        #   match[1] : the oxygen that bridges to phosphorus (this is our candidate connection point)
        #   match[2] : the phosphorus atom, etc.
        o_phospho = mol.GetAtomWithIdx(match[1])
        # Get neighbors of the oxygen excluding the phosphorus (match[2]) to get the glycerol carbon.
        glycerol_candidates = [nbr for nbr in o_phospho.GetNeighbors() if nbr.GetIdx() != match[2] and nbr.GetAtomicNum() == 6]
        if not glycerol_candidates:
            continue
        for sn3 in glycerol_candidates:
            # sn3 is expected to be the sn-3 carbon (a primary carbon, typically CH2).
            # Find a candidate sn-2: a carbon neighbor of sn3 (other than the original oxygen).
            sn3_neighbors = [atom for atom in sn3.GetNeighbors() if atom.GetIdx() != o_phospho.GetIdx() and atom.GetAtomicNum() == 6]
            for sn2 in sn3_neighbors:
                # Now, in glycerol, sn-2 (central carbon) should connect to both sn3 and sn1.
                sn2_neighbors = [atom for atom in sn2.GetNeighbors() if atom.GetIdx() != sn3.GetIdx() and atom.GetAtomicNum() == 6]
                if not sn2_neighbors:
                    continue
                for sn1 in sn2_neighbors:
                    # At this stage we have a candidate glycerol backbone: sn1 - sn2 - sn3.
                    # Now, check that each carbon bears the appropriate substituents.
                    
                    # 1. sn-2 should have an oxygen (other than those in the backbone) that leads to an acyl (ester) group.
                    acyl_found = False
                    for nbr in sn2.GetNeighbors():
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [sn1.GetIdx(), sn3.GetIdx()]:
                            if is_acyl_branch(nbr):
                                acyl_found = True
                                break
                    if not acyl_found:
                        continue
                        
                    # 2. sn-1 should have an oxygen (not linking to sn2) that leads to an alkyl (ether) chain.
                    alkyl_found = False
                    for nbr in sn1.GetNeighbors():
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != sn2.GetIdx():
                            if is_alkyl_branch(nbr):
                                alkyl_found = True
                                break
                    if not alkyl_found:
                        continue
                        
                    # 3. We already anchored sn-3 by its connection from phosphocholine.
                    # If the above conditions are met, accept the candidate.
                    return True, ("Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure "
                                  "with glycerol backbone (sn-1: alkyl, sn-2: acyl, sn-3: phosphocholine)")
    
    return False, "No glycerol backbone with required substituents (alkyl, acyl, and phosphocholine) found"

# Example usage (the examples below correspond to various structures in the class):
if __name__ == "__main__":
    test_smiles_list = [
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC(=O)[H])COCCCCCCCCCCCCCCCC",  # 1-hexadecyl-2-formyl-sn-glycero-3-phosphocholine
        "CCCCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCC\\C=C/C\\C=C/C\\C=C/CCCCC",  # 1-hexadecyl-2-[(8Z,11Z,14Z)-eicosatrienoyl]-sn-glycero-3-phosphocholine
        "P(OC[C@@H](COCCCCCCCCCCCCCCCC)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(=O)(OCC[N+](C)(C)C)[O-]"  # 1-O-hexadecyl-2-arachidonoyl-sn-glycero-3-phosphocholine
    ]
    for smi in test_smiles_list:
        result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")