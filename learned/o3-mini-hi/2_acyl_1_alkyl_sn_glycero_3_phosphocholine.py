"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
Definition: An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl
groups are located at positions 1 and 2 respectively.
Improved strategy:
  - Identify a phosphocholine headgroup: a phosphorus atom connected to at least one oxygen
    that is linked to a choline fragment (an N with a positive formal charge). Then among its
    oxygen neighbors, select one that connects to a carbon; this is taken as the tether (sn-3)
    to the glycerol backbone.
  - Traverse from the glycerol linkage (sn-3) to identify a contiguous three‐carbon chain.
    Identify the “middle” carbon (sn-2) and the “end” carbon (sn-1).
  - Confirm that sn-2 carries an acyl (ester) branch [an oxygen leading to a carbonyl] and
    sn-1 carries an alkyl (ether) branch (an oxygen branch that does not contain a carbonyl).
"""

from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based
    on its SMILES string. The method scans for a phosphocholine headgroup and then
    attempts to trace a glycerol backbone (sn-3 connecting to P, then sn-2 and sn-1),
    verifying that sn-2 bears an acyl (ester) branch and sn-1 bears an alkyl (ether) branch.

    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if a valid structure is found, False otherwise.
        str: Explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- HELPER FUNCTIONS ---
    def is_acyl_branch(o_atom):
        """
        Checks if an oxygen atom leads to an acyl branch.
        It does so by checking if the O connects to a carbon that has at least one
        double bond to oxygen (i.e. a carbonyl group).
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    def is_alkyl_branch(o_atom):
        """
        Checks if an oxygen atom leads to a pure alkyl branch.
        It verifies that an attached carbon is not part of a carbonyl group.
        """
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                has_carbonyl = False
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            has_carbonyl = True
                            break
                if not has_carbonyl:
                    return True
        return False

    # --- STEP 1: Identify the phosphocholine headgroup ---
    # Look for a phosphorus atom (atomic number 15) that has oxygen neighbors.
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    for ph in phos_atoms:
        oxy_neighbors = [nbr for nbr in ph.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            continue  # insufficient oxygen neighbors for phosphocholine

        # Identify the oxygen that leads to a choline fragment (i.e. attached to an N with positive charge)
        o_choline = None
        other_oxys = []
        for o in oxy_neighbors:
            attached_choline = False
            for nbr in o.GetNeighbors():
                if nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() > 0:
                    attached_choline = True
                    break
            if attached_choline:
                o_choline = o
            else:
                other_oxys.append(o)
        if o_choline is None or len(other_oxys) == 0:
            continue  # this phosphorus does not show a clear phosphocholine signature
        
        # Use one of the other oxygen neighbors as the linkage to the glycerol backbone.
        o_gly = other_oxys[0]
        # The oxygen must be attached to a carbon (the sn-3 carbon)
        gly_cands = [nbr for nbr in o_gly.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not gly_cands:
            continue  # no carbon found attached to the oxygen linking phosphate
        sn3 = gly_cands[0]
        
        # --- STEP 2: Trace the glycerol backbone along three connected carbons ---
        # In a glycerol backbone, sn-3 (bearing the P-linkage) connects to sn-2, which in turn connects to sn-1.
        # We look for a two-step carbon chain from sn3.
        sn3_neighbors = [nbr for nbr in sn3.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != o_gly.GetIdx()]
        if not sn3_neighbors:
            continue  # cannot find sn-2 candidate
        for sn2 in sn3_neighbors:
            sn2_neighbors = [nbr for nbr in sn2.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != sn3.GetIdx()]
            if not sn2_neighbors:
                continue  # cannot find sn-1 candidate
            for sn1 in sn2_neighbors:
                # Now check substituents on sn-2 and sn-1.
                # At sn-2: look for an oxygen (other than those bridging sn3 or sn1) that gives an acyl branch.
                sn2_oxygens = [nbr for nbr in sn2.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [sn3.GetIdx(), sn1.GetIdx()]]
                acyl_ok = False
                for ox in sn2_oxygens:
                    if is_acyl_branch(ox):
                        acyl_ok = True
                        break
                if not acyl_ok:
                    continue  # missing acyl branch at sn-2

                # At sn-1: look for an oxygen (other than that linking to sn2) that gives an alkyl branch.
                sn1_oxygens = [nbr for nbr in sn1.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != sn2.GetIdx()]
                alkyl_ok = False
                for ox in sn1_oxygens:
                    if is_alkyl_branch(ox):
                        alkyl_ok = True
                        break
                if not alkyl_ok:
                    continue  # missing alkyl branch at sn-1

                # If we reached here, we seem to have found a glycerol backbone linked at sn-3,
                # with sn-2 bearing an acyl (ester) branch and sn-1 bearing an alkyl (ether) branch.
                return True, ("Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine: "
                              "phosphocholine headgroup found with appropriate glycerol backbone "
                              "(sn-3 linked to phosphate; sn-2 has acyl branch; sn-1 has alkyl branch)")
    
    return False, "No valid phosphocholine-based glycerol backbone with required acyl (sn-2) and alkyl (sn-1) substituents found"

# Example usage; running some test SMILES (the provided examples) if run as a script.
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