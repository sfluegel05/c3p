"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid
Definition: Any terpenoid derived from a tetraterpene. The term includes compounds in which the C40 skeleton 
of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Heuristic criteria (updated):
  - Carbon count: between 20 and 120 carbons (to allow appended groups).
  - Extended conjugated polyene chain: defined as a series of at least 3 double bonds connected to each other.
    We ignore the built‐in conjugation flag and instead require that the bond is a double bond between carbon atoms.
  - Molecular weight: between 200 and 2000 Da.
Note: This heuristic approach may not capture every borderline structure.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string using heuristic criteria.
    
    Heuristics:
      1. The molecule must contain between 20 and 120 carbon atoms.
      2. The molecular weight must be between 200 and 2000 Da.
      3. The molecule must contain an extended conjugated polyene chain. In this code we search for 
         bonds that are double bonds between carbon atoms and build a connectivity graph between 
         them (two bonds are connected if they share an atom). We then compute the longest chain 
         (by number of bonds in a simple path) and require that it contains at least 3 bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a tetraterpenoid by heuristic, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbons)
    if c_count < 20 or c_count > 120:
        return False, f"Carbon count {c_count} outside acceptable range (20-120) for tetraterpenoids"
    
    # Check molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200 or mw > 2000:
        return False, f"Molecular weight {mw:.1f} out of range (200-2000 Da) for typical tetraterpenoids"
    
    # Identify double bonds between carbon atoms.
    # We ignore the 'conjugated' flag and simply check for C=C bonds.
    double_bond_indices = []
    bond_to_atoms = {}
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                bid = bond.GetIdx()
                double_bond_indices.append(bid)
                bond_to_atoms[bid] = (a1.GetIdx(), a2.GetIdx())
    
    if not double_bond_indices:
        return False, "No carbon–carbon double bonds found"
    
    # Build an adjacency graph among the double bonds.
    # Two double bonds are connected if they share an atom.
    graph = {bid: set() for bid in double_bond_indices}
    n = len(double_bond_indices)
    for i in range(n):
        bid1 = double_bond_indices[i]
        atoms1 = set(bond_to_atoms[bid1])
        for j in range(i+1, n):
            bid2 = double_bond_indices[j]
            atoms2 = set(bond_to_atoms[bid2])
            if atoms1.intersection(atoms2):
                graph[bid1].add(bid2)
                graph[bid2].add(bid1)
    
    # Use DFS to compute the longest simple path (number of bonds) in the bond graph.
    def dfs(current, visited):
        max_length = 0
        for neighbor in graph[current]:
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length
    
    longest_chain = 0
    for bid in double_bond_indices:
        chain_length = dfs(bid, {bid})
        if chain_length > longest_chain:
            longest_chain = chain_length
            
    # We require at least 3 connected double bonds.
    if longest_chain < 3:
        return False, "No sufficiently long conjugated polyene chain detected (need at least 3 connected C=C bonds)"
    
    return True, "Molecule meets carbon count, extended polyene chain, and molecular weight criteria typical for tetraterpenoids"

# Example usage:
if __name__ == "__main__":
    # Example SMILES string (feel free to substitute with one of the provided examples)
    example_smiles = "CO[C@@H]1[C@H](C)O[C@@H](O[C@@H](\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)C=C\\C=C(C)C=C\\C1)C(C)(C)O)[C@@H](OC)[C@@H]1O"
    result, reason = is_tetraterpenoid(example_smiles)
    print("Classified as tetraterpenoid?", result)
    print("Reason:", reason)