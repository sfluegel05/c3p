"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: Lipopolysaccharide

Definition:
    Lipopolysaccharides (LPS) are natural compounds (found in Gram-negative bacterial cell walls)
    typically consisting of a trisaccharide repeating unit (including two heptose units and 
    octulosonic acid) together with oligosaccharide side chains and multiple 3‐hydroxytetradecanoic 
    acid units. 

Heuristic criteria used in this classifier:
  1. The molecule should contain at least one sugar‐like ring (rings of size 5–7 with ≥1 oxygen).
     Additionally, it should contain at least one ring of size 7 or 8 (to mimic a heptose or 
     octulosonic acid unit).
  2. A carboxyl/carboxylic acid or ester motif is required; we search for a carbonyl bound 
     to an oxygen ([CX3](=O)[OX2]).
  3. A long aliphatic (nonaromatic carbon-only) chain is expected. Here we compute the longest 
     contiguous path of nonaromatic carbon atoms (allowing single or double bonds) and require it 
     to be at least 10 atoms long.
  4. Overall molecular weight is expected to be above ~800 Da.
  
Note:
  This is a heuristic method and will not capture the full complexity of LPS structures.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest contiguous chain consisting only of nonaromatic carbon atoms.
    Bonds considered are SINGLE or DOUBLE.
    """
    max_length = 0
    # List indices of carbon atoms that are not aromatic.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()]
    
    # Build adjacency between these carbons using only single or double bonds.
    adj = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
            continue
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        if idx1 in adj and idx2 in adj:
            adj[idx1].append(idx2)
            adj[idx2].append(idx1)
    
    # Depth-first search from each carbon atom.
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Criterion 1a: Look for at least one sugar-like ring (size 5–7 with ≥1 oxygen and the rest carbons)
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring_count = 0
    for ring in ring_info:
        if len(ring) not in [5, 6, 7]:
            continue
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        if symbols.count("O") >= 1 and symbols.count("C") >= len(ring) - 1:
            sugar_ring_count += 1
    if sugar_ring_count < 1:
        return False, f"Found only {sugar_ring_count} sugar-like ring(s); at least 1 is expected"

    # Criterion 1b: Check for a heptose (or octulosonic acid)–like ring (size 7 or 8 with ≥1 oxygen and mostly carbons)
    hepta_like_count = 0
    for ring in ring_info:
        if len(ring) in [7, 8]:
            symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            # Accept if at least one oxygen and rest predominantly carbon (i.e. number of carbons is at least len(ring)-1)
            if symbols.count("O") >= 1 and symbols.count("C") >= len(ring) - 1:
                hepta_like_count += 1
    if hepta_like_count < 1:
        return False, "No heptose/octulosonic acid–like ring (size 7 or 8) found, which is expected in LPS"

    # Criterion 2: Carboxyl/carboxylic acid or ester motif [CX3](=O)[OX2]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if acid_pattern is None:
        return False, "Error compiling acid SMARTS pattern"
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxyl/carboxylic acid (or ester) motif ([CX3](=O)[OX2]) found"

    # Criterion 3: Check longest contiguous nonaromatic carbon chain 
    longest_chain = longest_carbon_chain(mol)
    if longest_chain < 10:
        return False, f"Longest contiguous carbon chain is only {longest_chain} atoms long; at least 10 are expected"

    # Criterion 4: Molecular weight (lowered threshold to catch genuine LPS such as Mer-WF3010)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a lipopolysaccharide"
    
    return True, ("Structure contains sugar-like ring(s) (including a heptose/octulosonic acid–like ring), "
                  f"a carboxyl/ester motif, and a long aliphatic chain (longest chain = {longest_chain} carbons) "
                  f"with a molecular weight of {mol_wt:.1f} Da, consistent with a lipopolysaccharide.")

# Example usage (for testing):
if __name__ == "__main__":
    # Mer-WF3010 (previously false negative) is now tested.
    test_smiles = "O=C(OC1C(O)C2(OC(C1OC3OC(C(O)C(C3O)O)COC(=O)/C=C/C=C/C=C/C)CO)OCC=4C2=C(O)C=C(O)C4)/C=C/C=C/CC(O)/C(=C/C=C/CCC(CC)C)/C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print("Lipopolysaccharide classification:", result)
    print("Reason:", reason)