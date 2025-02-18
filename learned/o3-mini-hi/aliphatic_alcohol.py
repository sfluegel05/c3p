"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies:An aliphatic alcohol, defined as 'An alcohol derived from an aliphatic compound.'
Improvement rationale (as gleaned from previous outcomes):
 - Disqualify molecules that contain free carboxylic acid/carboxylate groups.
 - For each –OH group, require that it is attached to a saturated (sp³), non‐aromatic carbon.
 - From that candidate carbon (ignoring the –OH oxygen) perform a DFS search over contiguous 
   chain carbons (acyclic and non‐aromatic). In the DFS we allow a tolerance of up to 2 sp2 centers 
   (isolated) but disallow any carbonyl-like carbons (i.e. a carbon double bonded to oxygen).
 - If any adjoining chain (not counting the candidate carbon itself) has at least 6 carbons then
   the molecule is classified as an aliphatic alcohol.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol has at least one hydroxyl (-OH) group attached to a saturated (sp3),
    non‐aromatic carbon and that candidate carbon is connected to a contiguous, open-chain 
    aliphatic (acyclic) region (allowing up to 2 isolated sp2 centres) of at least 6 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Disqualify molecules that contain free carboxylic acid/carboxylate groups.
    acid_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[O;H1]"),  # carboxylic acid group
        Chem.MolFromSmarts("[CX3](=O)[O-]")     # carboxylate group
    ]
    for patt in acid_patterns:
        if mol.HasSubstructMatch(patt):
            return False, "Molecule contains carboxylic acid/carboxylate functionality"

    # --- Helper: determine if a carbon atom is allowed as part of the aliphatic chain.
    #    Returns "sp3" if the carbon is allowed as sp3, "sp2" if it is allowed as an isolated sp2, or False otherwise.
    def allowed_chain_atom(atom):
        if atom.GetAtomicNum() != 6:
            return False  # Only carbons allowed.
        # Disallow aromatic atoms and those in rings.
        if atom.GetIsAromatic() or atom.IsInRing():
            return False
        # Exclude carbons having a double bond to oxygen (carbonyl-like)
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return False
        # Allow sp3; allow sp2 only as isolated (tolerance will be enforced in DFS)
        hyb = atom.GetHybridization()
        if hyb == Chem.rdchem.HybridizationType.SP3:
            return "sp3"
        elif hyb == Chem.rdchem.HybridizationType.SP2:
            return "sp2"
        return False

    # --- DFS function to compute maximum contiguous chain length along allowed (acyclic, non‐aromatic) carbons.
    # We use backtracking with a visited set per branch and a tolerance counter for sp2 centres.
    def dfs(atom, visited, sp2_count):
        max_length = 0  # length counted from this starting atom.
        for nbr in atom.GetNeighbors():
            # Only consider carbon neighbors.
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue  # avoid cycles
            typ = allowed_chain_atom(nbr)
            if not typ:
                continue
            new_sp2 = sp2_count + (1 if typ == "sp2" else 0)
            if new_sp2 > 2:  # exceed allowed sp2 tolerance
                continue
            visited.add(nbr.GetIdx())
            current_length = 1 + dfs(nbr, visited, new_sp2)
            visited.remove(nbr.GetIdx())
            if current_length > max_length:
                max_length = current_length
        return max_length

    # --- Main search: iterate over oxygens to spot hydroxyl (-OH) groups.
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum() != 8:
            continue
        # Check that oxygen has at least one hydrogen (to signify -OH instead of, e.g., ether)
        if oxygen.GetTotalNumHs() < 1:
            continue

        # For each neighbor of the oxygen, check if it is a candidate carbon (bearing the –OH)
        for cand_carbon in oxygen.GetNeighbors():
            if cand_carbon.GetAtomicNum() != 6:
                continue
            # Must be saturated (sp3) and non-aromatic.
            if cand_carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or cand_carbon.GetIsAromatic():
                continue
            
            # For each neighbor of the candidate carbon (except the -OH oxygen), attempt to find a contiguous chain.
            for neigh in cand_carbon.GetNeighbors():
                if neigh.GetIdx() == oxygen.GetIdx():
                    continue  # Don't go back to the -OH oxygen.
                typ = allowed_chain_atom(neigh)
                if not typ:
                    continue
                # Start a DFS from this neighbor.
                visited = set([neigh.GetIdx()])
                chain_length = dfs(neigh, visited, 1 if typ=="sp2" else 0)
                # chain_length here is the count of chain carbons from the starting neighbor.
                if chain_length >= 6:
                    return True, ("Found -OH group attached to a sp3, non-aromatic carbon that is connected "
                                  "to an acyclic aliphatic chain (chain length = {}).".format(chain_length))
    return False, ("No qualifying -OH group found attached to a sp3, non-aromatic carbon with an adjoining "
                   "acyclic aliphatic chain of at least 6 carbons.")
    
# (Optional) Testing examples.
if __name__ == "__main__":
    test_cases = [
        # True positives from the provided outcomes:
        ("CCCCCCC(C)O", "octan-2-ol"),
        ("O=C(OC)/C=C/CC(O)CCCCCC(O)C", "Cladosporester A"),
        ("CCCCCC(O)CC(=O)CCc1ccc(O)c(O)c1", "1-(3,4-Dihydroxyphenyl)-5-hydroxy-3-decanone"),
        ("OCCCCCC/C=C/CCCC", "7E-dodecen-1-ol"),
        ("O[C@@H](CCC/C=C\\CC)C", "6Z-Nonen-2R-ol"),
        ("OCCCCCCCCCCC1=CC=CC=C1", "10-Phenyl-1-decanol"),
        ("CCCCCCCCCCCCCCCCCCCCCCO", "tricosan-1-ol"),
        # A few false positives/negatives reported before could be tested as well.
    ]
    
    for smi, name in test_cases:
        result, reason = is_aliphatic_alcohol(smi)
        print("SMILES:", smi)
        print("Name:", name)
        print("Classification:", result)
        print("Reason:", reason)
        print("-" * 60)