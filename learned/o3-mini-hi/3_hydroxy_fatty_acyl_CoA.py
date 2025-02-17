"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA resulting from the condensation of the thiol group 
of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.
Improved Strategy:
  1. Parse the molecule and add explicit hydrogens.
  2. Look for a thioester group (a carbonyl [C]=O directly attached to an S atom).
  3. For each such match, treat the carbonyl (C1) as the start of the acyl chain and then:
       - Find a candidate “α–carbon” (C2) attached to the carbonyl (excluding the carbonyl oxygen and the S).
       - From the α–carbon, find a neighboring carbon (candidate “β–carbon”, C3) that is sp3, is not in a ring, and bears at least one hydroxyl (–OH) group attached via a single bond.
       - Also perform a DFS along only aliphatic (acyclic) carbon–atoms starting from the α–carbon to ensure the fatty chain is of adequate length (at least two carbons beyond the carbonyl, i.e. α and β).
  4. Look for a diagnostic CoA fragment (we use a SMARTS that should match the phosphate–ribose–adenine part).
Note: This heuristic is not perfect. Real molecules or edge cases may deviate.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES.
    
    Requirements:
      1. Must contain an acyl thioester group (R-C(=O)-S-R').
      2. The fatty acyl group must be at least long enough so that, numbering the carbonyl as C1,
         there is an α–carbon (C2) and a β–carbon (C3) that bears an –OH group (the 3-hydroxy condition).
      3. Must contain a fragment diagnostic of Coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3-hydroxy fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for clarity.
    mol = Chem.AddHs(mol)
    
    # --- Step 1: Look for a thioester group ---
    # Our SMARTS pattern requires a carbonyl carbon (sp2) attached to a sulfur.
    thioester_smarts = "[CX3](=O)S"
    patt_thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(patt_thioester)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) functional group found"
    
    # Helper to check if an atom has a bonded -OH group.
    def has_OH(atom):
        # Look among neighbors for an oxygen attached via a single bond that in turn carries a hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Check if the oxygen has an explicit hydrogen (or implicit count > 0)
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Helper: Traverse an acyclic carbon chain (only sp3 carbons, not in any ring)
    # starting from a given atom and return the set of carbon atom indices reached.
    def traverse_chain(start_atom, exclude_idx):
        visited = set()
        stack = [start_atom]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            for nbr in atom.GetNeighbors():
                # Only follow carbon atoms (atomic number 6) that are not the excluded atom and not in a ring.
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != exclude_idx and not nbr.IsInRing():
                    if nbr.GetIdx() not in visited:
                        stack.append(nbr)
        return visited

    found_acyl = False
    # For each thioester match, look on the acyl (fatty) side.
    for match in thioester_matches:
        # In our pattern [CX3](=O)S, match[0] is the carbonyl carbon.
        carbonyl = mol.GetAtomWithIdx(match[0])
        # Look for neighbors that can be the beginning of the acyl chain.
        # Skip the carbonyl oxygen and the sulfur used in the thioester.
        acyl_candidates = []
        for nbr in carbonyl.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(carbonyl.GetIdx(), nbr.GetIdx())
            # Skip the double-bonded oxygen.
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                continue
            # Skip the sulfur (thioester linkage to CoA).
            if nbr.GetAtomicNum() == 16:
                continue
            if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()):
                acyl_candidates.append(nbr)
        if not acyl_candidates:
            continue
        
        # For each candidate α–carbon (which will be C2 in the acyl chain)
        for alpha in acyl_candidates:
            # From α–carbon, try to find a candidate “β–carbon” (C3) 
            # which is a neighbor of alpha (except going back to the carbonyl)
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == carbonyl.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()):
                    beta_candidates.append(nbr)
            if not beta_candidates:
                continue
            # For each beta candidate, check that it is an sp3 carbon and has a hydroxyl group.
            for beta in beta_candidates:
                if beta.GetHybridization() != rdchem.HybridizationType.SP3:
                    continue
                if not has_OH(beta):
                    continue
                # Also, ensure that the chain is at least as long as required.
                # Traverse the acyl chain starting from alpha (excluding the carbonyl)
                chain = traverse_chain(alpha, carbonyl.GetIdx())
                # Minimal fatty acyl chain (excluding the carbonyl) should have at least 2 carbons
                # (alpha and beta); here we use >=2.
                if len(chain) < 2:
                    continue
                # We accept this thioester as bearing a 3-hydroxy fatty acyl chain.
                found_acyl = True
                break
            if found_acyl:
                break
        if found_acyl:
            break

    if not found_acyl:
        return False, "No acyl chain found with an -OH at the β–position (3-hydroxy condition)"
    
    # --- Step 2: Look for a Coenzyme A fragment ---
    # We use a SMARTS that matches part of the CoA moiety (phosphate–ribose–adenine fragment).
    # Note: In some SMILES the phosphates may carry charges so we ignore the charge by using a generic pattern.
    coa_smarts = "COP(O)(=O)OP(O)(=O)OC"  
    patt_coa = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(patt_coa):
        return False, "Coenzyme A moiety not detected (using fragment pattern)"
    
    return True, "Contains a thioester group with a 3-hydroxy fatty acyl chain (OH on β–carbon) and a CoA moiety"


# If run as a script, perform some tests.
if __name__ == "__main__":
    # List (some examples only; many examples were in the specification)
    test_smiles = [
        # True positive examples:
        "CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # (R)-3-hydroxyhexanoyl-CoA
        "CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # 3-hydroxytetracosatetraenoyl-CoA
        # False positive example (should not be accepted):
        "C[C@@H]1O[C@@H](OCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)"
        "COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)N3C=NC4=C(N)N=CN=C34)[C@H](O)C[C@H]1O"  # Example known to be a false positive in previous attempt.
    ]
    
    for smi in test_smiles:
        result, reason = is_3_hydroxy_fatty_acyl_CoA(smi)
        print(result, "->", reason)