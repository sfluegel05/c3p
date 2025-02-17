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
       - The acyl chain side (ignoring the S side) should provide a candidate “α–carbon”.
       - Then look among the neighbors (other than the carbonyl) of the candidate α–carbon for a “β–carbon”
         that carries an –OH group.
       - Also, traverse the acyl chain (acyclic carbons only) to ensure that it has at least two carbons.
  4. Finally, require that a diagnostic Coenzyme A fragment is present.
Note: Because of unsaturation and chiral tagging, some filtering (such as checking sp3 hybridization)
      may be too strict; we have relaxed that check for the β–carbon.
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
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # Step 1: Identify the thioester group.
    # SMARTS: a carbon with three connections and a double-bonded oxygen, bonded to a sulfur.
    thioester_smarts = "[CX3](=O)S"
    patt_thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(patt_thioester)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) functional group found"
    
    # Helper function: check if an atom has at least one -OH neighbor.
    def has_OH(atom):
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # Ensure single bond (explicitly compare bond type)
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check if the oxygen has at least one hydrogen (explicit or implicit)
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Helper: Traverse an acyclic carbon chain starting from a given atom.
    # This ensures that the fatty acyl chain has at least two carbons (alpha and beta).
    def traverse_chain(start_atom, exclude_idx):
        visited = set()
        stack = [start_atom]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == exclude_idx:
                    continue
                if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                    if nbr.GetIdx() not in visited:
                        stack.append(nbr)
        return visited

    found_acyl = False
    # Loop over each thioester match.
    for match in thioester_matches:
        # In our SMARTS "[CX3](=O)S", match[0] is the carbonyl carbon and match[1] is the sulfur.
        carbonyl = mol.GetAtomWithIdx(match[0])
        
        # Find the neighbor on the acyl (fatty acid) side.
        acyl_candidates = []
        for nbr in carbonyl.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            # Skip the double-bonded oxygen (check using bond type).
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                continue
            # Skip the sulfur linked to CoA.
            if nbr.GetAtomicNum() == 16:
                continue
            # We expect an acyl chain carbon (not in ring).
            if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()):
                acyl_candidates.append(nbr)
        if not acyl_candidates:
            continue  # Try next thioester match
        
        # For each candidate α–carbon, check its neighbors for a β–carbon with an –OH.
        for alpha in acyl_candidates:
            # Get all neighbors (potential β–carbons) except going back to the carbonyl.
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == carbonyl.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()):
                    beta_candidates.append(nbr)
            if not beta_candidates:
                continue
            # Check each beta candidate for the presence of an –OH group.
            for beta in beta_candidates:
                # (Relaxed check: do not insist on sp3 so as to capture unsaturated cases.)
                if not has_OH(beta):
                    continue
                # Verify minimal chain length using a DFS from the α–carbon.
                chain_atoms = traverse_chain(alpha, carbonyl.GetIdx())
                if len(chain_atoms) < 2:
                    continue
                found_acyl = True
                break
            if found_acyl:
                break
        if found_acyl:
            break

    if not found_acyl:
        return False, "No acyl chain found with an -OH on the β–carbon (3-hydroxy condition)"
    
    # Step 2: Identify a diagnostic Coenzyme A fragment.
    # We use a SMARTS pattern matching part of the phosphate–ribose–adenine fragment.
    # This is a simplified pattern that should appear in the CoA moiety.
    coa_smarts = "COP(O)(=O)OP(O)(=O)OC"
    patt_coa = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(patt_coa):
        return False, "Coenzyme A moiety not detected (fragment pattern missing)"
    
    return True, "Contains a thioester group with a fatty acyl chain bearing an -OH on the β–carbon and a CoA moiety"


# If run as a script, perform some tests.
if __name__ == "__main__":
    # Test examples (a few from the list provided).
    test_smiles = [
        # True positives:
        "CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # (R)-3-hydroxyhexanoyl-CoA
        "CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # (3R,9Z,12Z,15Z,18Z)-3-hydroxytetracosatetraenoyl-CoA
        # (Some additional examples can be added from the specification)
        # False positive example:
        "C[C@@H]1O[C@@H](OCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)"
        "COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)N3C=NC4=C(N)N=CN=C34)[C@H](O)C[C@H]1O"
    ]
    
    for smi in test_smiles:
        result, reason = is_3_hydroxy_fatty_acyl_CoA(smi)
        print(result, "->", reason)