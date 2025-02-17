"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid
Definition: Any terpenoid derived from a tetraterpene. The term includes compounds in which the C40 skeleton 
of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Heuristic criteria (current version):
  1. Acceptable carbon count: between 20 and 120 carbon atoms.
  2. Molecular weight: between 200 and 2000 Da.
  3. Presence of an extended conjugated polyene chain: at least 3 consecutive conjugated double bonds.
     To detect this, we build a graph of atoms connected by conjugated double bonds and search for the longest
     simple path (measured as the number of bonds). This is to better cover stereochemical variants and avoid
     mis‐detection.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import BondType
import sys

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string using heuristic criteria.
    
    Heuristics:
      1. The molecule must contain between 20 and 120 carbon atoms.
      2. Its molecular weight must be between 200 and 2000 Da.
      3. The molecule must contain an extended conjugated polyene chain, defined here as a chain
         of at least three consecutive conjugated double bonds.
         We compute the longest such chain by building a graph of atoms linked by conjugated double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a tetraterpenoid by heuristic, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count < 20 or c_count > 120:
        return False, f"Carbon count {c_count} outside acceptable range (20-120) for tetraterpenoids"
    
    # Check molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200 or mw > 2000:
        return False, f"Molecular weight {mw:.1f} out of range (200-2000 Da) for typical tetraterpenoids"
    
    # Build a graph of atoms connected by conjugated double bonds.
    # key: atom index; value: list of tuples (neighbor_atom_index, bond index)
    double_bond_graph = {}
    for bond in mol.GetBonds():
        if bond.GetBondType() == BondType.DOUBLE and bond.GetIsConjugated():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            double_bond_graph.setdefault(a1, []).append((a2, bond.GetIdx()))
            double_bond_graph.setdefault(a2, []).append((a1, bond.GetIdx()))
            
    # If there are no conjugated double bonds at all then no polyene is present.
    if not double_bond_graph:
        return False, "No conjugated double bonds detected"
    
    # Depth-first search to find the longest chain (count as number of double bonds in the path)
    def dfs(current, prev, used_edges):
        max_length = len(used_edges)
        # For each neighbor reached through a conjugated double bond:
        for (nbr, bond_idx) in double_bond_graph.get(current, []):
            if bond_idx in used_edges:
                continue  # this bond was already used in the current path
            # In a chain the same edge should not be used twice.
            used_edges.add(bond_idx)
            # Continue DFS from the neighbor.
            current_length = dfs(nbr, current, used_edges)
            if current_length > max_length:
                max_length = current_length
            used_edges.remove(bond_idx)
        return max_length

    # Try starting DFS from every atom that is present in the double bond graph.
    longest_chain = 0
    for atom_idx in double_bond_graph.keys():
        for (nbr, bond_idx) in double_bond_graph[atom_idx]:
            used = set([bond_idx])
            chain_length = dfs(nbr, atom_idx, used)
            if chain_length > longest_chain:
                longest_chain = chain_length
            # Early exit if chain length requirement is met.
            if longest_chain >= 3:
                break
        if longest_chain >= 3:
            break

    if longest_chain < 3:
        return False, "No sufficiently long conjugated polyene chain detected (need at least 3 connected C=C bonds)"
    
    return True, "Molecule meets criteria for tetraterpenoid (carbon count, molecular weight, and polyene chain)"

# Example usage:
if __name__ == "__main__":
    # Replace 'example_smiles' with any SMILES string to test.
    example_smiles = "CO[C@@H]1[C@H](C)O[C@@H](O[C@@H](\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)C=C\\C=C(C)C=C\\C1)C(C)(C)O)[C@@H](OC)[C@@H]1O"
    result, reason = is_tetraterpenoid(example_smiles)
    print("Classified as tetraterpenoid?", result)
    print("Reason:", reason)