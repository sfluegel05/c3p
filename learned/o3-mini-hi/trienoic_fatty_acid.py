"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: A free (terminal) carboxylic acid with a long, acyclic (unbranched) acyl chain
that contains exactly three carbon–carbon double bonds. In addition, the three double bonds 
should occur in a methylene‐interrupted pattern (i.e. not directly conjugated).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    The criteria are:
      1. The molecule must be valid.
      2. It must contain exactly one free carboxylic acid group (SMARTS "[CX3](=O)[OX2H1]"),
         where the acid carbon is terminal (bonded to exactly one carbon).
      3. Starting from the acid carbon, the acyl chain is extracted by following the unique
         carbon neighbor (i.e. the chain must be linear and unbranched).
      4. The extracted (acyclic) chain must be sufficiently long (here at least 10 carbons)
         and contain exactly three C=C bonds.
      5. The double bonds must be “methylene‐interrupted” (i.e. no two double bonds occur
         on consecutive bonds along the fatty acyl chain).
      6. The molecular weight must be above a minimal threshold.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the trienoic fatty acid criteria, False otherwise.
        str: A short reason message.
    """
    # 1. Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Look for a free carboxylic acid group.
    # SMARTS: carbonyl carbon attached to [OX2H1]
    ca_smarts = "[CX3](=O)[OX2H1]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if len(ca_matches) == 0:
        return False, "Missing carboxylic acid functional group (COOH)"
    if len(ca_matches) > 1:
        return False, f"Found {len(ca_matches)} carboxylic acid groups; requires exactly 1"
    
    # In the SMARTS, the acid carbon is the first atom.
    acid_carbon_idx = ca_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    # 2a. Ensure the acid carbon is terminal: it must have exactly one carbon neighbor.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxylic acid group is not terminal (acid carbon must be bonded to exactly one carbon)"
    
    # 3. Extract the linear acyl chain.
    # Instead of a DFS that may pick a branched route, we follow the unique neighbor chain.
    # Typically, in fatty acids the chain is unbranched. (Branches lead to false hits.)
    chain_atoms = []  # will hold atom indices of the acyl chain (starting at the alpha carbon)
    double_bond_count = 0  # count bonds that are C=C (in the chain)
    double_bond_positions = []  # positions (in chain; bond between chain[i] and chain[i+1])
    
    # Start from the unique neighbor (first carbon of the chain)
    current_atom = carbon_neighbors[0]
    prev_atom = acid_carbon  # To avoid going back
    chain_atoms.append(current_atom.GetIdx())
    # Check the bond between acid group carbon and this carbon (normally not double)
    bond = mol.GetBondBetweenAtoms(prev_atom.GetIdx(), current_atom.GetIdx())
    # (We ignore any unsaturation here, assuming the acid carbon is part of the COOH group.)
    
    # Now iteratively follow the chain.
    pos = 0  # position in the chain list (0 means the bond from acid to first chain carbon is not counted)
    while True:
        # Get carbon neighbors of current_atom except the one we came from.
        nbrs = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()]
        # If more than one, the chain is branched.
        if len(nbrs) > 1:
            return False, f"Acyl chain is branched at atom index {current_atom.GetIdx()}; fatty acid chain must be linear"
        if len(nbrs) == 0:
            # End of chain reached.
            break
        # Otherwise, there is exactly one next carbon.
        next_atom = nbrs[0]
        # Check the bond between current_atom and next_atom
        bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), next_atom.GetIdx())
        if bond is None:
            return False, "Unexpected error: missing bond in acyl chain"
        # If it is a double bond, count it.
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bond_count += 1
            double_bond_positions.append(pos)
        # Append next_atom to the chain.
        chain_atoms.append(next_atom.GetIdx())
        # Update for next iteration.
        prev_atom = current_atom
        current_atom = next_atom
        pos += 1
    
    # 4. Require that the chain is long enough.
    if len(chain_atoms) < 10:
        return False, f"Acyl chain too short ({len(chain_atoms)} carbons); must be a long-chain fatty acid"
    
    # 5. The chain must be fully acyclic. (Given the linear follow above and no branching,
    #    we already reject ring participation via branching. But double-check: none of the atoms
    #    in the chain should be in a ring.)
    for idx in chain_atoms:
        if mol.GetAtomWithIdx(idx).IsInRing():
            return False, "Acyl chain contains rings; fatty acid chain must be acyclic"
    
    # 6. There must be exactly three carbon–carbon double bonds along the acyl chain.
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} C=C double bonds along the acyl chain; requires exactly 3 for a trienoic fatty acid"
    
    # 7. Verify that the double bonds occur in a methylene‐interrupted pattern.
    # That is, if d0, d1, d2 are the positions (in our chain, where bond between chain[i] and chain[i+1])
    # then we require that consecutive double bond positions differ by at least 2.
    if len(double_bond_positions) >= 2:
        for i in range(1, len(double_bond_positions)):
            if (double_bond_positions[i] - double_bond_positions[i-1]) < 2:
                return False, "Double bonds are conjugated; they must be separated by a methylene group"
    
    # 8. Overall molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for a fatty acid"
    
    return True, ("Contains one free, terminal carboxylic acid group and an acyclic, unbranched acyl chain "
                  "with exactly 3 C=C bonds in a methylene‐interrupted pattern typical of trienoic fatty acids")

# Example usage (for testing):
if __name__ == '__main__':
    # Test with one example: 10,12,14-octadecatrienoic acid.
    test_smiles = "OC(=O)CCCCCCCC/C=C/C=C/C=C/CCC"
    result, reason = is_trienoic_fatty_acid(test_smiles)
    print(result, reason)