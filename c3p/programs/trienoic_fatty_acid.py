"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: A polyunsaturated fatty acid that contains exactly three carbon–carbon double bonds,
with a terminal carboxylic acid (free acid, not esterified) and a long acyclic aliphatic chain.
This improved program extracts the main acyl chain starting at the free carboxyl carbon and 
counts double bonds strictly along that chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    The function verifies that:
      - The molecule is valid.
      - It contains exactly one free carboxylic acid group ([CX3](=O)[OX2H1]) whose acid carbon
        is terminal (bonded to only one carbon).
      - The main acyl chain (the longest contiguous chain of carbon atoms starting from the acid carbon)
        contains exactly three C=C bonds.
      - The acyl chain is long enough and lacks cyclic content.
      - The overall molecular weight is above a minimal threshold.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a trienoic fatty acid, False otherwise.
        str: A reason describing the classification result.
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Find free carboxylic acid group using SMARTS.
    # SMARTS for a free carboxylic acid (COOH) group.
    ca_smarts = "[CX3](=O)[OX2H1]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    
    if len(ca_matches) == 0:
        return False, "Missing carboxylic acid functional group (COOH)"
    if len(ca_matches) > 1:
        return False, f"Found {len(ca_matches)} carboxylic acid groups; requires exactly 1"

    # The match returns indices for the atoms in the pattern.
    # In our SMARTS "[CX3](=O)[OX2H1]", the first atom is the acid carbon.
    ca_indices = ca_matches[0]
    ca_carbon_idx = ca_indices[0]
    ca_carbon = mol.GetAtomWithIdx(ca_carbon_idx)
    
    # 2a. Ensure the acid carbon is terminal (bonded to exactly one carbon).
    carbon_neighbors = [nbr for nbr in ca_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon must be bonded to exactly one carbon)"
    
    # 3. Extract the main acyl chain by traversing only carbon atoms
    def get_longest_chain(mol, start_idx):
        """
        Returns the longest linear chain (as a list of atom indices) starting from start_idx.
        This DFS only follows carbon atoms and avoids backtracking.
        """
        def dfs(atom_idx, parent_idx):
            atom = mol.GetAtomWithIdx(atom_idx)
            paths = []
            extended = False
            # Only consider carbon neighbors (avoid the carboxyl oxygens etc).
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 6:
                    continue
                nbr_idx = nbr.GetIdx()
                if nbr_idx == parent_idx:
                    continue
                for sub_path in dfs(nbr_idx, atom_idx):
                    paths.append([atom_idx] + sub_path)
                    extended = True
            if not extended:
                return [[atom_idx]]
            return paths

        all_paths = dfs(start_idx, -1)
        # Choose the path with maximum length.
        longest = max(all_paths, key=len)
        return longest

    chain = get_longest_chain(mol, ca_carbon_idx)
    # For a fatty acid, we expect the chain to be fairly long; here we require at least 10 carbons.
    if len(chain) < 10:
        return False, f"Acyl chain too short ({len(chain)} carbons); must be a long-chain fatty acid"
    
    # 4. Ensure that none of the atoms in the chain are part of any ring.
    for idx in chain:
        if mol.GetAtomWithIdx(idx).IsInRing():
            return False, "Acyl chain contains rings; fatty acid chain must be acyclic"
    
    # 5. Count the number of carbon–carbon double bonds along the main chain.
    double_bond_count = 0
    for i in range(len(chain)-1):
        bond = mol.GetBondBetweenAtoms(chain[i], chain[i+1])
        if bond is None:
            # if for some reason the path is broken, skip.
            continue
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} C=C double bonds along the acyl chain; requires exactly 3 for a trienoic fatty acid"
    
    # 6. Check overall molecular weight to avoid very small molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a fatty acid"
    
    # Passed all the tests, so it qualifies as a trienoic fatty acid.
    return True, ("Contains a free, terminal carboxylic acid group and a long acyclic acyl chain with exactly 3 "
                  "carbon–carbon double bonds typical of trienoic fatty acids")

# Example usage:
if __name__ == '__main__':
    # Test with one known trienoic fatty acid: 10,12,14-octadecatrienoic acid.
    test_smiles = "OC(=O)CCCCCCCC/C=C/C=C/C=C/CCC"
    result, reason = is_trienoic_fatty_acid(test_smiles)
    print(result, reason)