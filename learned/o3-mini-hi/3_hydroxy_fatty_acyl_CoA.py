"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
Definition: A hydroxy fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.
Improved strategy:
  1. Parse the molecule and add explicit hydrogens.
  2. Look for a thioester group – a carbonyl (C(=O)) attached to an S atom.
  3. From the carbonyl atom, identify the acyl chain (skip the S and the C=O oxygen) and require:
      • the chain to have (when numbering from the carbonyl as C1) an identifiable α carbon and then a β carbon,
      • the β carbon (C3 overall) carries a hydroxyl (–OH) that is bound by a single bond.
      • also, we require the acyl chain to contain at least 3 carbon atoms (a minimal “fatty” chain).
  4. Look for a fragment diagnostic of Coenzyme A. Here we use a pattern that
     matches the phosphate–ribose–adenine part “COP(O)(=O)OP(O)(=O)OC” that is typically present.
     
Note: This heuristic is not perfect. Real molecules may deviate from these patterns.
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    Requirements:
      1. It must contain an acyl thioester bond (i.e. R-C(=O)-S-R').
      2. The fatty acyl group (R) must be at least three carbons long such that,
         counting the carbonyl carbon as position 1, the β carbon (C3) carries an -OH.
      3. It must contain a fragment diagnostic of Coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-hydroxy fatty acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # --- Step 1: Identify a thioester group ---
    # We look for a pattern: a carbonyl (sp2 carbon with =O) directly attached to a sulfur.
    thioester_smarts = "[CX3](=O)S"
    patt_thioester = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(patt_thioester)
    if not thioester_matches:
        return False, "No thioester (acyl-CoA) functional group found"
    
    # Helper: check that an atom has a bonded -OH group.
    def has_OH(atom):
        # Look at neighbor atoms for an O that is singly bonded and itself bonded to a hydrogen.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                # Look for explicit hydrogen on that oxygen.
                for onbr in nbr.GetNeighbors():
                    if onbr.GetAtomicNum() == 1:
                        return True
        return False

    found_hydroxy = False
    acyl_chain_found = False
    # For each thioester match, inspect the “acyl” side.
    for match in thioester_matches:
        # In our pattern [CX3](=O)S,
        # match[0] is the carbonyl carbon.
        carbonyl_atom = mol.GetAtomWithIdx(match[0])
        # From the carbonyl, find a neighbor that is carbon (the beginning of the acyl chain)
        acyl_candidates = []
        for nbr in carbonyl_atom.GetNeighbors():
            # Exclude the oxygen in the carbonyl:
            bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                continue
            # Exclude the sulfur (the thioester linkage to CoA)
            if nbr.GetAtomicNum() == 16:
                continue
            if nbr.GetAtomicNum() == 6:  # carbon on acyl chain
                acyl_candidates.append(nbr)
        if not acyl_candidates:
            continue  # no suitable acyl chain found for this thioester

        # For each candidate starting carbon (this will be our α-carbon, i.e. C2 when carbonyl is C1)
        for alpha in acyl_candidates:
            # To reduce false positives, require that the candidate is not in a ring
            if alpha.IsInRing():
                continue
            # From the α-carbon, try to find a neighboring carbon (other than the carbonyl)
            beta_candidates = []
            for nbr in alpha.GetNeighbors():
                if nbr.GetIdx() == carbonyl_atom.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6:
                    beta_candidates.append(nbr)
            if not beta_candidates:
                continue
            # For each possible β-carbon candidate, require that:
            # 1) The β carbon is not in a ring (i.e. part of a free acyl chain)
            # 2) It carries an -OH group (as a -CH(OH)- would have).
            for beta in beta_candidates:
                if beta.IsInRing():
                    continue
                if has_OH(beta):
                    # Now, do a simple check that the acyl branch is of at least minimal length.
                    # Perform a DFS starting from the α-carbon (exclude the carbonyl) following only carbon atoms.
                    visited = set()
                    stack = [alpha]
                    chain_carbons = set()
                    while stack:
                        atom = stack.pop()
                        if atom.GetIdx() in visited:
                            continue
                        visited.add(atom.GetIdx())
                        if atom.GetAtomicNum() == 6:
                            chain_carbons.add(atom.GetIdx())
                        for nbr in atom.GetNeighbors():
                            # Only follow carbons that are not the carbonyl (which is not part of the chain)
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_atom.GetIdx():
                                if nbr.GetIdx() not in visited:
                                    stack.append(nbr)
                    # Require at least 3 carbon atoms in the acyl chain (i.e. C1 from carbonyl not counted,
                    # then alpha and beta at minimum).
                    if len(chain_carbons) < 2:
                        continue
                    found_hydroxy = True
                    acyl_chain_found = True
                    break
            if found_hydroxy:
                break
        if found_hydroxy:
            break

    if not found_hydroxy:
        return False, "No acyl chain with an OH at the 3-position was found"
    
    # --- Step 2: Identify a Coenzyme A fragment ---
    # We use a more explicit CoA fragment SMARTS. In many acyl-CoA molecules the CoA part
    # contains the fragment "COP(O)(=O)OP(O)(=O)OC" (part of the adenosine 3’,5’-bisphosphate diphosphate).
    coa_smarts = "COP(O)(=O)OP(O)(=O)OC"
    patt_coa = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(patt_coa):
        return False, "Coenzyme A moiety not detected (using refined pattern)"
    
    return True, "Contains a thioester group with a 3-hydroxy fatty acyl chain (OH on β-carbon) and a CoA moiety"


# Testing examples (the valid examples have been given in the specifications)
if __name__ == "__main__":
    test_smiles = [
        "CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # (R)-3-hydroxyhexanoyl-CoA, valid
        "C[C@@H]1O[C@@H](OCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)"
        "COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)N3C=NC4=C(N)N=CN=C34)[C@H](O)C[C@H]1O"  # false positive in previous attempt
    ]
    for smi in test_smiles:
        result, reason = is_3_hydroxy_fatty_acyl_CoA(smi)
        print(result, "->", reason)