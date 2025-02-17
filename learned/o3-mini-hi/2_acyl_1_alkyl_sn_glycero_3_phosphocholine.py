"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
Definition: An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl
groups are located at positions 1 and 2 respectively.
"""

from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    The algorithm is as follows:
      1. Parse the molecule.
      2. Verify that a phosphocholine headgroup is present using a SMARTS filter.
      3. Identify a glycerol backbone substructure defined as a contiguous three-carbon chain, 
         ideally CH2–CH–CH2.
      4. For the three carbons (sn-1, sn-2, sn-3) check that a non-backbone oxygen is attached:
             • sn-1: an oxygen whose branch is an alkyl (ether) chain (ideally long).
             • sn-2: an oxygen whose branch is part of an ester (acyl) group (i.e. attached to C(=O)...).
             • sn-3: an oxygen whose branch is the phosphocholine group (O–P with P attached to a trimethylammonium).
      5. If such a backbone is found with the correct substituents, return True.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is in the class, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check that a phosphocholine headgroup fragment is present.
    # This SMARTS targets an oxygen attached to a phosphorus which in turn is linked to a positively charged nitrogen.
    phospho_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phospho_frag = Chem.MolFromSmarts(phospho_smarts)
    if not mol.HasSubstructMatch(phospho_frag):
        return False, "Phosphocholine headgroup not found"
    
    # Helper function to check if an oxygen leads to a phosphocholine branch.
    def is_phospho_branch(o_atom):
        # Look for a phosphorus neighbor connected to a nitrogen with positive formal charge.
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 15:  # phosphorus
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 7 and subnbr.GetFormalCharge() == 1:
                        return True
        return False

    # Helper function to check if an oxygen leads to an acyl (ester) branch
    def is_acyl_branch(o_atom):
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                # look for a double bond from this carbon to an oxygen (C=O)
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # Helper function to check if an oxygen leads to an alkyl branch.
    # Here we traverse a short fragment and count the number of carbon atoms.
    def is_alkyl_branch(o_atom):
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # Start a simple DFS search from this neighbor (avoid going back to oxygen o_atom)
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
                    # Limit search depth to avoid going everywhere.
                    if c_count > 5:
                        return True
                    for nn in current.GetNeighbors():
                        # Do not go back to the oxygen branch or into heteroatoms.
                        if nn.GetAtomicNum() == 6 and nn.GetIdx() not in visited:
                            stack.append(nn)
        return False

    # Look for a glycerol backbone.
    # We define a glycerol backbone as three connected carbons with pattern CH2-CH-CH2.
    # The SMARTS below matches a chain of three carbons.
    glycerol_smarts = "[CH2]-[CH]-[CH2]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    
    if not matches:
        return False, "No three-carbon backbone (potential glycerol) found"
    
    # For each potential glycerol backbone, test the three positions.
    for match in matches:
        # match is a tuple of atom indices corresponding to the three carbons.
        c1 = mol.GetAtomWithIdx(match[0])
        c2 = mol.GetAtomWithIdx(match[1])
        c3 = mol.GetAtomWithIdx(match[2])
        
        # For each carbon in the backbone, we need to find a non-backbone oxygen neighbor.
        # We record which branch matches which expected function.
        branch_found = {"alkyl": False, "acyl": False, "phospho": False}
        
        # Helper: given a carbon and expected branch function, find if one non-backbone oxygen qualifies.
        def check_branch(carbon, branch_type):
            for nbr in carbon.GetNeighbors():
                # Skip if neighbor is part of the backbone.
                if nbr.GetIdx() in match:
                    continue
                if nbr.GetAtomicNum() == 8:  # oxygen
                    if branch_type == "phospho" and is_phospho_branch(nbr):
                        return True
                    if branch_type == "acyl" and is_acyl_branch(nbr):
                        return True
                    if branch_type == "alkyl" and is_alkyl_branch(nbr):
                        return True
            return False

        # According to nomenclature for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine:
        # sn-1 carbon (c1) should carry the alkyl branch.
        # sn-2 carbon (c2) should carry the acyl branch.
        # sn-3 carbon (c3) should carry the phosphocholine branch.
        if not check_branch(c1, "alkyl"):
            continue  # try next backbone
        if not check_branch(c2, "acyl"):
            continue
        if not check_branch(c3, "phospho"):
            continue
        
        return True, ("Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure "
                      "(glycerol backbone with defined alkyl, acyl, and phosphocholine branches)")
    
    return False, "No glycerol backbone with the required substituents (alkyl, acyl, and phosphocholine) found"

# Example usage (you can remove or comment out these lines when deploying the module):
if __name__ == "__main__":
    test_smiles_list = [
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC(=O)[H])COCCCCCCCCCCCCCCCC",  # 1-hexadecyl-2-formyl-sn-glycero-3-phosphocholine
        "CCCCCCCCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCC\\C=C/C\\C=C/C\\C=C/CCCCC"  # example with eicosatrienoyl branch
    ]
    for smi in test_smiles_list:
        result, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")