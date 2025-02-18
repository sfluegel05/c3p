"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds,
with a free (terminal) carboxylic acid group and a long, acyclic acyl chain.
This program extracts the longest contiguous carbon chain starting from the acid carbon,
requires that the acid carbon is terminal, and then counts the C=C bonds along that chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    Checks that:
    - The molecule is valid.
    - It contains exactly one free carboxylic acid group ([CX3](=O)[OX2H1]),
      and that the acid carbon is terminal (bonded to exactly one other carbon).
    - The longest acyl chain (extracted starting from the acid carbon using DFS with
      a visited set) is sufficiently long, acyclic, and contains exactly three C=C bonds.
    - The overall molecular weight meets a minimal threshold.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule fits the trienoic fatty acid definition, False otherwise.
        str: A reason explaining the classification result.
    """
    # 1. Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Locate free carboxylic acid group: uses SMARTS "[CX3](=O)[OX2H1]".
    ca_smarts = "[CX3](=O)[OX2H1]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    
    if len(ca_matches) == 0:
        return False, "Missing carboxylic acid functional group (COOH)"
    if len(ca_matches) > 1:
        return False, f"Found {len(ca_matches)} carboxylic acid groups; requires exactly 1"

    # In the carboxylic acid SMARTS, the first atom is the acid carbon.
    ca_indices = ca_matches[0]
    acid_carbon_idx = ca_indices[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # 2a. Ensure the acid carbon is terminal (should have exactly one carbon neighbor).
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon must be bonded to exactly one carbon)"

    # 3. Extract the longest contiguous acyl chain following only carbon atoms (no backtracking).
    def get_longest_chain(mol, start_idx):
        """
        Uses depth-first search (DFS) with a visited list to find the longest linear chain 
        (list of atom indices) that begins at the start_idx. Only carbon atoms are followed.
        """
        longest_chain = []
        def dfs(current_idx, path):
            nonlocal longest_chain
            path.append(current_idx)
            # Update longest_chain if this path is longer.
            if len(path) > len(longest_chain):
                longest_chain = list(path)
            # Look through neighbors that are carbons, but do not revisit atoms on current path.
            for nbr in mol.GetAtomWithIdx(current_idx).GetNeighbors():
                if nbr.GetAtomicNum() != 6:
                    continue
                nbr_idx = nbr.GetIdx()
                if nbr_idx in path:
                    continue  # avoid cycles
                dfs(nbr_idx, path)
            path.pop()  # backtrack
        
        dfs(start_idx, [])
        return longest_chain

    chain = get_longest_chain(mol, acid_carbon_idx)
    # For a fatty acid, we expect a long chain (here requiring at least 10 carbon atoms).
    if len(chain) < 10:
        return False, f"Acyl chain too short ({len(chain)} carbons); must be a long-chain fatty acid"
    
    # 4. Verify that no atoms in the acyl chain are part of any ring.
    for idx in chain:
        if mol.GetAtomWithIdx(idx).IsInRing():
            return False, "Acyl chain contains rings; fatty acid chain must be acyclic"
    
    # 5. Count the number of carbon–carbon double bonds along the acyl chain.
    double_bond_count = 0
    for i in range(len(chain)-1):
        bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
        if bond is None:
            continue  # safety check; every sequential pair must be bonded
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} C=C double bonds along the acyl chain; requires exactly 3 for a trienoic fatty acid"
    
    # 6. Check overall molecular weight - a minimal threshold to avoid very small molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a fatty acid"
    
    return True, ("Contains one free, terminal carboxylic acid group and a long acyclic acyl chain "
                  "with exactly 3 carbon–carbon double bonds typical of trienoic fatty acids")

# Example usage (for testing):
if __name__ == '__main__':
    # Example: 10,12,14-octadecatrienoic acid.
    test_smiles = "OC(=O)CCCCCCCC/C=C/C=C/C=C/CCC"
    result, reason = is_trienoic_fatty_acid(test_smiles)
    print(result, reason)