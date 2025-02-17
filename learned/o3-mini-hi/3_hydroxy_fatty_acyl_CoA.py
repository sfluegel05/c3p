"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA resulting from the condensation of the thiol group 
of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.
Improved strategy:
  1. Parse the molecule and add explicit hydrogens.
  2. Look for a thioester group (a C(=O)S fragment).
  3. For each thioester match, treat the carbonyl carbon (C1) and then examine its acyl (fatty acid) neighbor:
       - Among neighbors of the carbonyl, the acyl chain side (exclude the S and carbonyl oxygen) is taken as the α–carbon (C2).
       - Then, among the neighbors of this candidate α–carbon (except the carbonyl), at least one must be a β–carbon (C3) bearing an –OH group.
       - Also, ensure that the acyl chain indeed has at least two carbons.
  4. Finally, require that a diagnostic Coenzyme A fragment is present.
"""

from rdkit import Chem
from rdkit.Chem import rdchem, AllChem, rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES.
    
    Requirements:
      1. Must contain a thioester group, i.e. a carbonyl (C(=O)) directly bonded to a sulfur.
      2. The fatty acyl chain (the side attached to the carbonyl, not the CoA portion) 
         must be long enough such that, labeling the carbonyl as C1, there is an α–carbon (C2)
         and a β–carbon (C3) that bears an –OH group.
      3. The molecule must contain a fragment diagnostic of Coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the 3-hydroxy fatty acyl-CoA heuristic, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES and check validity
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are available
    mol = Chem.AddHs(mol)
    
    # Step 1: Identify the thioester group using SMARTS.
    # Pattern: a tetrahedral or trigonal planar carbon with a double bonded O and a single-bond to S.
    thioester_smarts = "[CX3](=O)S"
    patt_thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(patt_thioester)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) functional group found"
    
    # Helper function: Checks if given atom bears an -OH (i.e. an oxygen attached by a single bond that also has at least one H)
    def beta_has_OH(atom):
        # Scan neighbors of atom for oxygen with at least one hydrogen and a single bond linkage.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() == rdchem.BondType.SINGLE and nbr.GetTotalNumHs() > 0:
                    return True
        return False
    
    found_acyl = False
    # Loop through each thioester match
    for match in thioester_matches:
        # In the SMARTS "[CX3](=O)S", match[0] is the carbonyl carbon; match[1] is the sulfur.
        carbonyl = mol.GetAtomWithIdx(match[0])
        
        # Identify neighbor atoms of carbonyl, excluding the oxygen from the C=O.
        acyl_candidates = []
        for nbr in carbonyl.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            # Skip the double-bonded oxygen
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                continue
            # Skip sulfur (the CoA leaving part)
            if nbr.GetAtomicNum() == 16:
                continue
            # The candidate acyl chain atom should be carbon and not in a ring.
            if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                acyl_candidates.append(nbr)
        if not acyl_candidates:
            # This thioester match does not offer an acyl (fatty acid) side.
            continue
        
        # For each candidate α-carbon, check for a β-carbon with an -OH group.
        for alpha in acyl_candidates:
            # For a proper fatty acyl chain we need at least one neighboring carbon (other than carbonyl)
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == carbonyl.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                    beta_candidates.append(nbr)
            if not beta_candidates:
                continue  # Not enough chain length
            
            # Check each candidate beta for an -OH substituent.
            for beta in beta_candidates:
                if beta_has_OH(beta):
                    found_acyl = True
                    break
            if found_acyl:
                break
        if found_acyl:
            break

    if not found_acyl:
        return False, "No acyl chain found with an -OH on the β–carbon (3-hydroxy condition)"
    
    # Step 2: Detect a diagnostic Coenzyme A fragment.
    # This simplified pattern searches for the phosphate-ribose portion.
    coa_smarts = "COP(O)(=O)OP(O)(=O)OC"
    patt_coa = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(patt_coa):
        return False, "Coenzyme A moiety not detected (fragment pattern missing)"
    
    return True, "Contains a thioester group with an acyl chain bearing an -OH on the β–carbon and a CoA moiety"

# If run as a script, perform some tests.
if __name__ == "__main__":
    # Test examples from the provided list
    test_smiles = [
        # (R)-3-hydroxyhexanoyl-CoA
        "CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        # (3R,9Z,12Z,15Z,18Z)-3-hydroxytetracosatetraenoyl-CoA
        "CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        # (R)-carnitinyl-CoA betaine
        "CC(C)(COP([O-])(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C[C@@H](O)C[N+](C)(C)C"
    ]
    
    for smi in test_smiles:
        result, reason = is_3_hydroxy_fatty_acyl_CoA(smi)
        print(result, "->", reason)