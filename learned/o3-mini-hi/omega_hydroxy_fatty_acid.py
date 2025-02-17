"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: Omega-hydroxy fatty acid
Definition: Any member of the class of naturally-occurring straight-chain fatty acids n carbon atoms long
with a carboxyl group (COOH) at position 1 and a hydroxyl (OH) at the omega (last) carbon.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an omega-hydroxy fatty acid.
    This requires a linear (unbranched) carbon chain with a terminal carboxylic acid group and an
    -OH group on the opposite terminal (omega) carbon.
  
    Args:
        smiles (str): SMILES string for the molecule.
  
    Returns:
        bool: True if the molecule meets the criteria for an omega-hydroxy fatty acid, False otherwise.
        str: Explanation for the determination.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to ensure all -OH groups are represented explicitly.
    mol = Chem.AddHs(mol)
    
    # Search for the carboxylic acid group.
    # Modify the SMARTS pattern to better capture a COOH group.
    # "[CX3](=O)[OX2H]" describes a tetra-coordinate carbon (sp2, in C=O)
    # attached to an -OH group (where the oxygen shows 2 connectivity [OX2] and has one hydrogen).
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (-COOH) group found in the molecule"

    # Assume the first match is the carboxyl group of interest.
    # The first index in the match tuple corresponds to the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Traverse all carbon atoms connected to the acid carbon to identify the fatty acid chain.
    carbon_indices = set()
    to_visit = [acid_carbon_idx]
    
    while to_visit:
        current = to_visit.pop()
        if current in carbon_indices:
            continue
        carbon_indices.add(current)
        atom = mol.GetAtomWithIdx(current)
        # Only consider neighbors that are carbons.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                if nbr.GetIdx() not in carbon_indices:
                    to_visit.append(nbr.GetIdx())
    
    if len(carbon_indices) < 2:
        return False, "The carbon chain is too short to be a fatty acid"
    
    # Check if the carbon chain is linear (unbranched).
    # In a linear chain, exactly two carbons (the termini) should have only one neighbor in the chain.
    carbon_degree = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        # Count only neighbors that belong to the carbon subgraph.
        deg = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbon_indices)
        carbon_degree[idx] = deg

    terminals = [idx for idx, deg in carbon_degree.items() if deg == 1]
    if len(terminals) != 2:
        return False, f"Carbon chain is not linear; expected 2 terminal carbons but found {len(terminals)}"
    
    # One terminal carbon should be the carboxyl carbon.
    if acid_carbon_idx not in terminals:
        return False, "Carboxylic acid group is not located at a terminal carbon"
    
    # The other terminal carbon is the omega candidate.
    omega_candidate_idx = [idx for idx in terminals if idx != acid_carbon_idx][0]
    omega_candidate = mol.GetAtomWithIdx(omega_candidate_idx)
    
    # Check for a hydroxyl (-OH) substituent on the omega candidate.
    has_omega_oh = False
    for nbr in omega_candidate.GetNeighbors():
        if nbr.GetAtomicNum() == 8:  # Check oxygen atoms
            # Count attached hydrogens (should have at least one)
            num_h = sum(1 for subnbr in nbr.GetNeighbors() if subnbr.GetAtomicNum() == 1)
            if num_h >= 1:
                has_omega_oh = True
                break

    if not has_omega_oh:
        return False, "Terminal (omega) carbon lacks a hydroxyl (-OH) substituent"

    return True, "Molecule is a linear fatty acid with a COOH group at C1 and an -OH at the omega position"


# Example test (remove or comment out these test calls as needed):
if __name__ == '__main__':
    # Example: (2E)-18-hydroxyoctadec-2-enoic acid
    test_smiles = "OCCCCCCCCCCCCCCC\\C=C\\C(O)=O"
    result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
    print(result, reason)