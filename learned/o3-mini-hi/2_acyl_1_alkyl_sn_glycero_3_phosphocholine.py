"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
#!/usr/bin/env python
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
Definition: An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl
groups are located at positions 1 and 2 respectively.
Improved strategy:
  - Identify the phosphocholine headgroup by scanning phosphorus atoms that are bound to three oxygens,
    one of which connects to a positively charged (choline) fragment.
  - Use the non-choline oxygen as the anchor to the glycerol backbone.
  - Trace the glycerol backbone (sn-3 → sn-2 → sn-1) and check that:
      sn-2 carries an acyl (ester) branch (i.e. an O–C(=O)– fragment)
      sn-1 carries an alkyl (ether) branch (i.e. an O–alkyl fragment without a carbonyl)
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
    
    # --- HELPER FUNCTIONS ---
    def is_acyl_branch(o_atom):
        """
        Checks if an oxygen atom leads to an acyl ester branch (should encounter a carbon bearing a C=O).
        """
        # Look at the neighbor that is carbon and check for a carbonyl bond.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                # Check bonds on this carbon for a double bond to oxygen
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    def is_alkyl_branch(o_atom, min_carbons=1):
        """
        Checks if an oxygen atom leads to a pure alkyl branch (does not contain a carbonyl).
        Uses a DFS to count the number of carbons and detects any C=O.
        """
        visited = set()
        stack = []
        # Start with carbons attached to the oxygen.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                stack.append(nbr)
        c_count = 0
        while stack:
            current = stack.pop()
            if current.GetIdx() in visited:
                continue
            visited.add(current.GetIdx())
            # If a double bond to oxygen (carbonyl) is found, then not a pure alkyl branch.
            for b in current.GetBonds():
                if b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other = b.GetOtherAtom(current)
                    if other.GetAtomicNum() == 8:
                        return False
            # Count the carbon.
            c_count += 1
            if c_count >= min_carbons:
                # We have enough carbon units for a minimal alkyl chain.
                return True
            # Continue exploring the branch.
            for nn in current.GetNeighbors():
                if nn.GetAtomicNum() == 6 and nn.GetIdx() not in visited:
                    stack.append(nn)
        return False

    # --- STEP 1: Identify the phosphocholine headgroup by scanning for a phosphorus (P) atom ---
    # Phosphocholine phosphorus should be bound to three oxygens: one involved in a double bond, one
    # perhaps carrying a negative charge, and two single-bonded oxygens. One of these oxygens should lead
    # to a choline fragment (look for an N with a positive formal charge).
    phos_candidates = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]  # phosphorus atoms
    for ph in phos_candidates:
        # Check if phosphorus has 3 oxygen neighbors (the fourth bond might be an additional attachment).
        neighbors = ph.GetNeighbors()
        oxy_neighbors = [at for at in neighbors if at.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            continue  # Too few oxygen neighbors to be phosphocholine
        # Verify that at least one bond to oxygen is a double bond (the phosphoryl oxygen).
        has_dbl = False
        for bond in ph.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtom(ph).GetAtomicNum() == 8:
                has_dbl = True
                break
        if not has_dbl:
            continue
        
        # Among oxygen neighbors, identify which one attaches to a choline fragment.
        possible_choline = []
        possible_gly = []
        for o in oxy_neighbors:
            # Look for an attached nitrogen with a positive charge.
            attached_choline = False
            for nbr in o.GetNeighbors():
                if nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() > 0:
                    attached_choline = True
                    break
            if attached_choline:
                possible_choline.append(o)
            else:
                possible_gly.append(o)
        if len(possible_choline) == 0 or len(possible_gly) == 0:
            continue  # Cannot identify both choline and glycerol linkage
        
        # Assume one of the non-choline oxygens (possible_gly) is the linkage to the glycerol backbone (sn-3).
        o_gly = possible_gly[0]
        
        # --- STEP 2: Trace to the glycerol backbone
        # From o_gly, the connecting atom must be a carbon (sn-3).
        sn3_candidates = [nbr for nbr in o_gly.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not sn3_candidates:
            continue  # no carbon attached to the glycerol oxygen
        for sn3 in sn3_candidates:
            # From sn3, the glycerol backbone should proceed to sn-2.
            # We require a neighboring carbon (other than the oxygen o_gly) as candidate for sn-2.
            sn3_neighbors = [nbr for nbr in sn3.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != o_gly.GetIdx()]
            if not sn3_neighbors:
                continue
            for sn2 in sn3_neighbors:
                # From sn2, expect at least one other carbon (sn-1) that is not sn3.
                sn2_neighbors = [nbr for nbr in sn2.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != sn3.GetIdx()]
                if not sn2_neighbors:
                    continue
                for sn1 in sn2_neighbors:
                    # --- Check substituents on the backbone positions ---
                    # (a) At sn-2: there should be an oxygen (excluding those linking to sn3 or sn1)
                    #     that leads to an acyl branch displaying a carbonyl (C=O) (i.e. ester).
                    acyl_found = False
                    for nbr in sn2.GetNeighbors():
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [sn3.GetIdx(), sn1.GetIdx()]:
                            if is_acyl_branch(nbr):
                                acyl_found = True
                                break
                    if not acyl_found:
                        continue  # missing acyl branch at sn-2
                    
                    # (b) At sn-1: there should be an oxygen (excluding the one linking back to sn2)
                    #     that leads to an alkyl (ether) substituent.
                    alkyl_found = False
                    for nbr in sn1.GetNeighbors():
                        if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != sn2.GetIdx():
                            if is_alkyl_branch(nbr):
                                alkyl_found = True
                                break
                    if not alkyl_found:
                        continue  # missing alkyl branch at sn-1
                    
                    # (c) Optionally, ensure the backbone connectivity is not excessive.
                    # For glycerol, sn-3 should be primarily connected to sn-2 and the linking oxygen.
                    sn3_carbon_neighbors = [a for a in sn3.GetNeighbors() if a.GetAtomicNum() == 6]
                    if len(sn3_carbon_neighbors) < 1:  # too many extraneous bonds
                        continue
                    
                    # If we have a valid glycerol skeleton with the proper branches, then we assume a match.
                    return True, ("Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine: glycerol backbone found "
                                  "(sn-3 linked to phosphocholine, sn-2 carries acyl branch, sn-1 carries alkyl branch)")
    
    return False, "No valid phosphocholine-based glycerol backbone with required acyl (sn-2) and alkyl (sn-1) substituents found"

# Example usage: testing a few SMILES strings from the provided list.
if __name__ == "__main__":
    test_smiles_list = [
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC(=O)[H])COCCCCCCCCCCCCCCCC",  # 1-hexadecyl-2-formyl-sn-glycero-3-phosphocholine
        "CCCCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCC\\C=C/C\\C=C/C\\C=C/CCCCC",  # 1-hexadecyl-2-[(8Z,11Z,14Z)-eicosatrienoyl]-sn-glycero-3-phosphocholine
        "CCCCCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCC\\C=C/C\\C=C/CCCCC",  # 1-octadecyl-2-[(13Z,16Z)-docosadienoyl]-sn-glycero-3-phosphocholine
        "P(OC[C@@H](COCCCCCCCCCCCCCCCC)OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(=O)(OCC[N+](C)(C)C)[O-]",  # 1-O-hexadecyl-2-arachidonoyl-sn-glycero-3-phosphocholine
        "CCCCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCC",  # 1-hexadecyl-2-butanoyl-sn-glycero-3-phosphocholine
    ]
    
    for smi in test_smiles_list:
        result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")