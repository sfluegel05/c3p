"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: Lipopolysaccharide

Definition:
  Lipopolysaccharides (LPS) are natural compounds found in Gram‐negative bacterial cell walls.
  They are built around a trisaccharide repeating unit – typically consisting of two heptose units
  and one octulosonic acid – together with oligosaccharide side chains and multiple 3‐hydroxytetradecanoic acid units.

Heuristic criteria used in this classifier:
  1. The molecule must contain at least three isolated sugar‐like rings. Here we define “sugar‐like”
     as a ring of size 5–7 that contains at least one oxygen; we further require that each ring is isolated,
     meaning that none of its atoms belongs to another ring.
  2. At least three carboxyl/carboxylic acid or ester motifs ([CX3](=O)[OX2]) must be found.
  3. There must be at least one long contiguous chain of nonaromatic carbon atoms (only single or double bonds)
     of length at least 14 (consistent with the 3‐hydroxytetradecanoic acid unit).
  4. Overall molecular weight should be at least ~800 Da.
  5. The ratio of oxygen to carbon atoms should be reasonably high (here ≥0.3) as expected in a saccharide-rich molecule.

Note:
  This is a heuristic method and will not capture the full complexity of LPS structures.
  
The approach was modified in response to previous false positive and false negative outcomes.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def longest_carbon_chain(mol):
    """
    Computes the length of the longest contiguous chain consisting only of nonaromatic carbon atoms.
    Only bonds of type SINGLE or DOUBLE are considered.
    """
    max_length = 0
    # List indices of carbon atoms that are not aromatic.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()]
    
    # Build an adjacency list between these carbons using only single or double bonds.
    adj = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
            continue
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        if idx1 in adj and idx2 in adj:
            adj[idx1].append(idx2)
            adj[idx2].append(idx1)
    
    def dfs(current, visited):
        length = 1  # count current atom
        local_max = length
        for neighbor in adj[current]:
            if neighbor not in visited:
                new_visited = visited | {neighbor}
                candidate = dfs(neighbor, new_visited)
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
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    # Count how many rings each atom is in.
    atom_ring_count = {i: 0 for i in range(mol.GetNumAtoms())}
    for ring in ring_info:
        for idx in ring:
            atom_ring_count[idx] += 1

    # Criterion 1: Count isolated sugar‐like rings.
    sugar_rings = 0
    for ring in ring_info:
        if len(ring) not in [5, 6, 7]:
            continue
        # Get atom symbols in the ring.
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        # A sugar-like ring should have at least one oxygen.
        if symbols.count("O") < 1:
            continue
        # Check that the ring is "isolated" (none of its atoms is in >1 ring).
        if all(atom_ring_count[idx] == 1 for idx in ring):
            sugar_rings += 1
    if sugar_rings < 3:
        return False, f"Found only {sugar_rings} isolated sugar-like ring(s); at least 3 are expected"
    
    # Criterion 2: Look for carboxyl/carboxylic acid or ester motif.
    # The SMARTS pattern [CX3](=O)[OX2] covers both acids and esters.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if acid_pattern is None:
        return False, "Error compiling carboxyl/ester SMARTS pattern"
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) < 3:
        return False, f"Found only {len(acid_matches)} carboxyl/ester motif(s); at least 3 are expected"
    
    # Criterion 3: Check for a long, contiguous nonaromatic carbon chain.
    chain_length = longest_carbon_chain(mol)
    if chain_length < 14:
        return False, f"Longest contiguous nonaromatic carbon chain is only {chain_length} atoms; at least 14 are expected"
    
    # Criterion 4: Check for overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); expected at least 800 Da"
    
    # Criterion 5: Check the oxygen to carbon ratio.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count == 0 or (o_count / c_count) < 0.3:
        return False, f"Oxygen/carbon ratio too low ({o_count}/{c_count}); expected at least 0.3"
    
    return True, (f"Structure contains {sugar_rings} isolated sugar-like ring(s), {len(acid_matches)} carboxyl/ester motif(s), "
                  f"a longest nonaromatic carbon chain of {chain_length} atoms, a molecular weight of {mol_wt:.1f} Da, and an O/C ratio of {o_count/c_count:.2f}, "
                  "which is consistent with a lipopolysaccharide.")

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one example SMILES (e.g. Mer-WF3010).
    test_smiles = "O=C(OC1C(O)C2(OC(C1OC3OC(C(O)C(C3O)O)COC(=O)/C=C/C=C/C=C/C)CO)OCC=4C2=C(O)C=C(O)C4)/C=C/C=C/CC(O)/C(=C/C=C/CCC(CC)C)/C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print("Lipopolysaccharide classification:", result)
    print("Reason:", reason)