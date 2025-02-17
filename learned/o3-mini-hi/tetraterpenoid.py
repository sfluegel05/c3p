"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid
Definition: Any terpenoid derived from a tetraterpene. The term includes compounds in which the C40 skeleton 
of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Heuristic criteria:
  - Relaxed carbon count: typically between 20 and 120 carbons (allowing modification via appended groups).
  - Presence of an extended conjugated polyene chain. Rather than a fixed SMARTS, we compute the longest 
    contiguous chain of conjugated C=C bonds (we require at least three such bonds).
  - Optional: check on molecular weight typical (e.g. between 200 and 1500 Da).
Note: This is a heuristic approach and may not work for every borderline structure.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    A tetraterpenoid is derived from a tetraterpene (nominally with a C40 skeleton) that may have been modified.
    Many such molecules (e.g., carotenoids) display long conjugated polyene chains.
    
    This function uses several heuristic criteria:
      1. The molecule should have a reasonable number of carbon atoms (here, 20-120).
      2. The molecule should contain an extended conjugated polyene chain. We compute the longest chain of
         consecutively connected, conjugated double bonds and require at least 3 such bonds.
      3. The molecular weight should be in a typical range for tetraterpenoids (here, between 200 and 1500 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a tetraterpenoid, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    # Allow a relaxed range for tetraterpenoid derivatives.
    if c_count < 20 or c_count > 120:
        return False, f"Carbon count {c_count} outside acceptable range (20-120) for tetraterpenoids"
    
    # Check molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200 or mw > 1500:
        return False, f"Molecular weight {mw:.1f} out of range (200-1500 Da) for typical tetraterpenoids"
    
    # Build a list of bonds that are double and conjugated.
    double_bond_ids = []
    bond_to_atoms = {}  # map bond index -> tuple(begin_atom_idx, end_atom_idx)
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetIsConjugated():
            idx = bond.GetIdx()
            double_bond_ids.append(idx)
            beg = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            bond_to_atoms[idx] = (beg, end)
    
    # Build an adjacency list among these bonds: two bonds are connected if they share an atom.
    graph = {bid: set() for bid in double_bond_ids}
    for i in range(len(double_bond_ids)):
        for j in range(i+1, len(double_bond_ids)):
            bid1 = double_bond_ids[i]
            bid2 = double_bond_ids[j]
            atoms1 = set(bond_to_atoms[bid1])
            atoms2 = set(bond_to_atoms[bid2])
            if atoms1.intersection(atoms2):
                graph[bid1].add(bid2)
                graph[bid2].add(bid1)
    
    # Compute the length of the longest simple path (in terms of number of bonds) in this graph.
    # Because these chains are usually short, a DFS from each node is acceptable.
    def dfs(current, visited):
        max_length = 0
        for neighbor in graph[current]:
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length

    longest_chain = 0
    for bid in double_bond_ids:
        chain_len = dfs(bid, {bid})
        if chain_len > longest_chain:
            longest_chain = chain_len

    # We require at least 3 consecutive conjugated double bonds.
    if longest_chain < 3:
        return False, "No sufficiently long conjugated polyene chain detected (need at least 3 connected double bonds)"
    
    return True, "Molecule meets carbon count, extended conjugated polyene chain, and molecular weight criteria typical for tetraterpenoids"

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES (feel free to change this example)
    example_smiles = "CO[C@@H]1[C@H](C)O[C@@H](O[C@@H](\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)C=C\\C=C(C)C=C\\C1)C(C)(C)O)[C@@H](OC)[C@@H]1O"
    result, reason = is_tetraterpenoid(example_smiles)
    print("Classified as tetraterpenoid?", result)
    print("Reason:", reason)