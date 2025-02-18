"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile 
Definition: Any nitrile derived from an aliphatic compound.
In this version we require that (i) the nitrile carbon is not aromatic;
(ii) its sole substituent (the R branch, i.e. the neighbor that is not the nitrile nitrogen)
    is a carbon atom; (iii) if that branch is “simple” (i.e. consists only of carbons and R is primary)
    and is longer than a short alkyl chain (more than 3 carbons) then we reject it;
(iv) we also disqualify a branch whose immediate neighbors indicate a strongly polar (functionalized)
group, such as a carbonyl.
This is a heuristic and may miss some complex cases.
"""

from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    
    Procedure:
      1. Parse the SMILES and find nitrile groups matching [C;!a;X2]#[N;X1].
      2. For each nitrile, take the nitrile carbon (must be non‐aromatic) and then identify its
         unique substituent R (the neighbor that is not the nitrile nitrogen). R must be a carbon.
      3. Evaluate the R substituent branch:
         a. If R is “simple” – meaning that R has only one heavy‐atom neighbor (i.e. is primary) 
            and if the entire branch (found by a breadth‐first search into non‐aromatic parts) 
            consists solely of carbon atoms and contains >3 carbons – consider it too trivial and reject.
         b. Also check that among the atoms directly attached to R (other than the nitrile carbon)
            no polar “carbonyl” (i.e. an oxygen doubly bonded to something) is present.
      4. If at least one nitrile qualifies, return True plus a message.
         Otherwise, return False with a brief explanation.
         
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an aliphatic nitrile, False otherwise.
        str: Explanation text.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a nitrile group: a non‐aromatic carbon with two connections triple-bonded to a nitrogen.
    nitrile_pattern = Chem.MolFromSmarts("[C;!a;X2]#[N;X1]")
    if nitrile_pattern is None:
        return False, "Error creating nitrile SMARTS pattern"
    
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile group found in the molecule"
    
    # For each nitrile group, examine its substituent branch.
    for match in matches:
        nitrile_c = mol.GetAtomWithIdx(match[0])
        nitrile_n = mol.GetAtomWithIdx(match[1])
        
        # (A) Ensure the nitrile carbon is non‐aromatic ("aliphatic part").
        if nitrile_c.GetIsAromatic():
            continue
        
        # (B) Identify the unique substituent R (neighbor that is not the nitrile nitrogen)
        neighbors = [nbr for nbr in nitrile_c.GetNeighbors() if nbr.GetIdx() != nitrile_n.GetIdx()]
        if not neighbors:
            continue  # unusual nitrile with no other substituent
        R = neighbors[0]
        # Require that R is a carbon atom.
        if R.GetSymbol() != "C":
            continue
        
        # (C) Check that R has at least one hydrogen (to indicate a saturated, non-quaternary site)
        if R.GetTotalNumHs() < 1:
            continue

        # (D) Check immediate neighbors of R (other than nitrile_c) for interfering polar groups.
        bad_attachment = False
        for nbr in R.GetNeighbors():
            if nbr.GetIdx() == nitrile_c.GetIdx():
                continue
            # If neighbor is oxygen then check if it is double-bonded (as in a carbonyl).
            if nbr.GetSymbol() == "O":
                for bond in mol.GetBonds():
                    # Look for a bond involving this oxygen that is a double bond
                    if (bond.GetBeginAtom().GetIdx() == nbr.GetIdx() or bond.GetEndAtom().GetIdx() == nbr.GetIdx()) \
                       and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        bad_attachment = True
                        break
            if bad_attachment:
                break
        if bad_attachment:
            continue
        
        # (E) Traverse the branch starting at R (not going back to nitrile_c)
        # Here we collect all atoms reachable from R (excluding nitrile_c) that are not aromatic.
        branch_atoms = set()
        queue = [R]
        while queue:
            atom = queue.pop(0)
            if atom.GetIdx() in branch_atoms:
                continue
            branch_atoms.add(atom.GetIdx())
            # Only traverse nonaromatic atoms; we do not want to pass into aromatic portions.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() == nitrile_c.GetIdx():
                    continue
                if not nbr.GetIsAromatic():
                    queue.append(nbr)
                    
        # Count how many carbon atoms are in the branch.
        branch_carbon_count = sum(1 for idx in branch_atoms if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
        
        # (F) Determine if the branch is “simple” (only carbons) and if so, whether it is too long.
        # We consider R to be "primary" if besides the nitrile connection, it has only one other heavy-atom neighbor.
        heavy_neighbors_R = [nbr for nbr in R.GetNeighbors() if nbr.GetIdx() != nitrile_c.GetIdx() and nbr.GetAtomicNum() > 1]
        is_primary = (len(heavy_neighbors_R) == 1)
        
        branch_is_simple = is_primary and all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in branch_atoms)
        # If the branch is a plain alkyl chain and is longer than 3 carbons, we decide it is too trivial.
        if branch_is_simple and branch_carbon_count > 3:
            continue
        
        # If we make it here, this nitrile group passes our tests.
        return True, "Contains a nitrile group attached to an exclusively aliphatic substituent branch"
    
    return False, "Nitrile group(s) found but none have a qualifying aliphatic substituent branch"

# Example usages (uncomment to test):
# print(is_aliphatic_nitrile("N#CCC#N"))          # malononitrile -> expected True.
# print(is_aliphatic_nitrile("OC(=O)CNCC#N"))      # N-(cyanomethyl)glycine -> expected True.
# print(is_aliphatic_nitrile("CCOC(=O)C(C)C#N"))   # ethyl 2-cyanopropionate -> expected False.
# print(is_aliphatic_nitrile("CC#N"))              # acetonitrile -> expected True.
# print(is_aliphatic_nitrile("OCC#N"))             # glycolonitrile -> expected False.