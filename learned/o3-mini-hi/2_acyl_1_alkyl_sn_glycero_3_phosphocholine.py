"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
Definition: An alkyl,acyl-sn-glycero-3-phosphocholine in which unspecified alkyl and acyl
groups are located at positions 1 and 2 respectively.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    The algorithm is as follows:
      1. Parse the molecule.
      2. Verify that a phosphocholine head group is present using a SMARTS filter.
      3. Search for a glycerol-like carbon (a carbon with three oxygen neighbours) that can act as the
         sn-glycerol “central” carbon.
      4. For that candidate, classify each oxygen substituent:
           • One oxygen must be attached to a phosphorus which in turn is bonded to a positively charged nitrogen 
             (i.e. a phosphocholine branch).
           • One oxygen must lead to an acyl ester: that is, the oxygen is connected to a carbon that bears a 
             carbonyl group (C(=O)).
           • The remaining oxygen is assumed to be the alkyl (ether) branch.
      5. If exactly one such glycerol center is found with the expected three branches we return True.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the class, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check for the phosphocholine headgroup.
    # This SMARTS targets an oxygen attached to a phosphorus that is bound to an ethylene linked trimethylammonium.
    phospho_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phospho_frag = Chem.MolFromSmarts(phospho_smarts)
    if not mol.HasSubstructMatch(phospho_frag):
        return False, "Phosphocholine headgroup not found"
    
    # Next, we look for the glycerol-like backbone.
    # We expect a glycerol central carbon to be an sp3 carbon with three oxygen neighbours.
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            # Count neighboring oxygens.
            o_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if len(o_neighbors) == 3:
                candidates.append((atom, o_neighbors))
                
    if not candidates:
        return False, "No glycerol-like carbon (with three oxygen substituents) found"

    # Try to classify one candidate as the glycerol central carbon with the proper 1-,2-,3-substitutions.
    # We expect exactly:
    #   • one oxygen branch attached via phosphorus -> phosphocholine (position 3)
    #   • one oxygen branch attached to a carbonyl group (ester) -> acyl group (position 2)
    #   • one remaining oxygen branch that is ether-linked -> alkyl group (position 1)
    for cent_atom, o_list in candidates:
        count_phosphate = 0
        count_acyl = 0
        count_alkyl = 0
        for o in o_list:
            assigned = False
            # Check if this oxygen is connected to a phosphorus (for phosphocholine)
            for nbr in o.GetNeighbors():
                if nbr.GetAtomicNum() == 15 and nbr.GetIdx() != cent_atom.GetIdx():
                    # Further check that phosphorus is attached to a positively charged nitrogen.
                    found_n_plus = False
                    for pnbr in nbr.GetNeighbors():
                        if pnbr.GetAtomicNum() == 7 and pnbr.GetFormalCharge() == 1:
                            found_n_plus = True
                            break
                    if found_n_plus:
                        count_phosphate += 1
                        assigned = True
                        break  # go to next oxygen branch
            if assigned:
                continue
            
            # Check if this oxygen leads to an acyl chain (ester bond).
            # This is the case if the oxygen is bonded to a carbon (not the glycerol central) that is carbonylated.
            for nbr in o.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != cent_atom.GetIdx():
                    # Look for a double bond from this carbon to an oxygen.
                    for bond in nbr.GetBonds():
                        # For ester, we expect a C(=O) group.
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            other = bond.GetOtherAtom(nbr)
                            if other.GetAtomicNum() == 8:
                                count_acyl += 1
                                assigned = True
                                break
                    if assigned:
                        break
            if assigned:
                continue

            # If not phosphate or acyl, we assume this oxygen leads to the alkyl branch.
            count_alkyl += 1

        if count_phosphate == 1 and count_acyl == 1 and count_alkyl == 1:
            return True, "Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure (glycerol backbone with phosphocholine, acyl, and alkyl branches)"

    return False, "No glycerol carbon with the required three substituents (phosphocholine, acyl, and alkyl) was identified"


# Example usage (you can remove or comment out these lines when deploying the module):
if __name__ == "__main__":
    # Test with one of the examples: 1-hexadecyl-2-formyl-sn-glycero-3-phosphocholine
    test_smiles = "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC(=O)[H])COCCCCCCCCCCCCCCCC"
    is_class, reason = is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(test_smiles)
    print(is_class, reason)