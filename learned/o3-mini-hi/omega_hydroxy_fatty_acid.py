"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: Omega-hydroxy fatty acid
Definition: A naturally-occurring straight-chain fatty acid with n carbon atoms that has a carboxyl group (COOH) at position 1 and a hydroxyl group (OH) at the omega (last) carbon.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an omega-hydroxy fatty acid.
    This requires a linear (unbranched) carbon chain with a terminal carboxylic acid group and an
    -OH group on the opposite terminal carbon.
  
    Args:
        smiles (str): SMILES string for the molecule.
  
    Returns:
        bool: True if the molecule meets the criteria for an omega-hydroxy fatty acid, False otherwise.
        str: Explanation for the determination.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that -OH groups are explicit
    mol = Chem.AddHs(mol)
    
    # Search for the carboxylic acid group.
    # The SMARTS "[CX3](=O)[OX1H]" matches a carbonyl carbon (sp2) bonded to an -OH.
    acid_smarts = "[CX3](=O)[OX1H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (-COOH) group found in the molecule"
    
    # For our purposes, we consider the first match.
    # In the match tuple, the first index corresponds to the carbon of the acid group.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Next, we wish to analyze the carbon chain that makes up the fatty acid.
    # We build the subgraph of all carbon atoms connected to the acid carbon.
    carbon_indices = set()
    to_visit = [acid_carbon_idx]
    
    while to_visit:
        current = to_visit.pop()
        if current in carbon_indices:
            continue
        carbon_indices.add(current)
        atom = mol.GetAtomWithIdx(current)
        # only traverse neighbors that are carbons (atomic number 6)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                if nbr.GetIdx() not in carbon_indices:
                    to_visit.append(nbr.GetIdx())
    
    if len(carbon_indices) < 2:
        return False, "The carbon chain is too short to be a fatty acid"
    
    # Now, check if the carbon chain is linear (i.e. forms an unbranched chain).
    # We count, for each carbon in the subgraph, how many neighbors (within the subgraph) it has.
    carbon_degree = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        # count only neighbors that are in our carbon chain
        deg = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbon_indices)
        carbon_degree[idx] = deg

    # In a linear (straight-chain) structure, there will be exactly two terminal carbons (degree 1)
    terminals = [idx for idx, deg in carbon_degree.items() if deg == 1]
    if len(terminals) != 2:
        return False, f"Carbon chain is not linear; expected 2 terminal carbons but found {len(terminals)}"
    
    # One terminal carbon should be the acid group carbon.
    if acid_carbon_idx not in terminals:
        return False, "Carboxylic acid group is not located at a terminal carbon"
    
    # The other terminal carbon is the omega (terminal) carbon candidate.
    omega_candidate_idx = [idx for idx in terminals if idx != acid_carbon_idx][0]
    omega_candidate = mol.GetAtomWithIdx(omega_candidate_idx)
    
    # Now, check if the omega candidate is substituted with a hydroxyl (-OH).
    # Look at its neighbors for an oxygen that carries at least one hydrogen.
    has_omega_oh = False
    for nbr in omega_candidate.GetNeighbors():
        if nbr.GetAtomicNum() == 8:
            # Check if oxygen has at least one attached hydrogen
            # After AddHs, explicit hydrogens are present; we can count them.
            num_h = sum(1 for subnbr in nbr.GetNeighbors() if subnbr.GetAtomicNum() == 1)
            if num_h >= 1:
                has_omega_oh = True
                break

    if not has_omega_oh:
        return False, "Terminal (omega) carbon lacks a hydroxyl (-OH) substituent"
    
    return True, "Molecule is a linear fatty acid with a COOH group at C1 and an -OH at the omega position"
    
# Example use (you can remove or comment out these test calls if integrating elsewhere):
if __name__ == '__main__':
    test_smiles = "OCCCCCCCCCCCCCCC\\C=C\\C(O)=O"  # (2E)-18-hydroxyoctadec-2-enoic acid
    result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
    print(result, reason)