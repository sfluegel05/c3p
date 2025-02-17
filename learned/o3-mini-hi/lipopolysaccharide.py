"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: Lipopolysaccharide

Definition:
    Lipopolysaccharides are natural compounds consisting of a trisaccharide repeating unit 
    (two heptose units and octulosonic acid) with oligosaccharide side chains and 
    3-hydroxytetradecanoic acid units. They are a major constituent of the cell walls
    of Gram-negative bacteria.

Heuristic criteria used in this classifier:
  1. At least one sugar‐like ring is detected. Here we look at rings of size 5–7 that contain at least one oxygen and all other atoms are mostly carbon.
  2. Presence of a carboxyl/carboxylic acid or ester motif is required. We look for a carbonyl bound to an oxygen ([CX3](=O)[OX2]).
  3. Presence of a long aliphatic chain is expected. In the previous version a SMARTS for 8 contiguous “C” atoms was used but that misses chains with unsaturations. Here we use a custom function (longest_carbon_chain) that does a depth‐first search along nonaromatic carbon atoms (allowing both single and double bonds) to estimate the longest “chain.”
  4. Overall molecular weight is expected to be >1000 Da.

Note:
  This is a heuristic approach. Lipopolysaccharides are very complex and structure‐based identification remains challenging.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest contiguous chain consisting only of carbon atoms
    (ignoring aromatic atoms) in the molecule. Bonds of type SINGLE or DOUBLE are considered.
    This function does a DFS from every carbon atom.
    """
    max_length = 0
    # Get list of indices for nonaromatic carbon atoms.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()]
    
    # Build an adjacency dictionary for nonaromatic carbon atoms linked by single or double bonds.
    adj = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        # Only consider bonds that are SINGLE or DOUBLE
        if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        if idx1 in adj and idx2 in adj:
            adj[idx1].append(idx2)
            adj[idx2].append(idx1)
    
    # DFS routine – starting at a carbon atom, record path length.
    def dfs(current, visited):
        length = 1  # count the current atom
        local_max = length
        for nbr in adj[current]:
            if nbr not in visited:
                new_visited = visited | {nbr}
                candidate = dfs(nbr, new_visited)
                if length + candidate > local_max:
                    local_max = length + candidate
        return local_max

    for idx in carbon_idxs:
        chain_length = dfs(idx, {idx})
        if chain_length > max_length:
            max_length = chain_length
    return max_length

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a lipopolysaccharide, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Check for at least one sugar-like ring.
    # We inspect rings of size 5-7, requiring at least one oxygen and that all others are carbons.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring_count = 0
    for ring in ring_info:
        if len(ring) not in [5, 6, 7]:
            continue
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        oxygen_count = symbols.count("O")
        carbon_count = symbols.count("C")
        # A typical sugar ring has one oxygen and the rest are carbons.
        if oxygen_count >= 1 and carbon_count >= (len(ring) - 1):
            sugar_ring_count += 1
    if sugar_ring_count < 1:
        return False, f"Found only {sugar_ring_count} sugar-like ring(s); at least 1 is expected"
    
    # Criterion 2: Check for a carboxyl/carboxylic acid or ester motif.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if acid_pattern is None:
        return False, "Error compiling acid SMARTS pattern"
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxyl/carboxylic acid (or ester) motif ([CX3](=O)[OX2]) found"
    
    # Criterion 3: Check for a long aliphatic chain.
    # We now compute the longest path of contiguous nonaromatic carbon atoms (allowing single or double bonds).
    longest_chain = longest_carbon_chain(mol)
    if longest_chain < 8:
        return False, f"Longest contiguous carbon chain is only {longest_chain} atoms long; at least 8 are expected"
    
    # Criterion 4: Check the overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a lipopolysaccharide"
    
    return True, ("Structure contains sugar-like ring(s), a carboxyl/ester motif, and a long aliphatic chain "
                  f"(longest chain = {longest_chain} carbons) with a molecular weight of {mol_wt:.1f} Da, "
                  "consistent with a lipopolysaccharide.")


# Example usage (for testing):
if __name__ == "__main__":
    # Test with one of the provided SMILES (Mer-WF3010) as an example.
    test_smiles = "O=C(OC1C(O)C2(OC(C1OC3OC(C(O)C(C3O)O)COC(=O)/C=C/C=C/C=C/C)CO)OCC=4C2=C(O)C=C(O)C4)/C=C/C=C/CC(O)/C(=C/C=C/CCC(CC)C)/C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print("Lipopolysaccharide classification:", result)
    print("Reason:", reason)