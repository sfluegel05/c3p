"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-hydroxy fatty acid

Definition:
A 2-hydroxy fatty acid is an acyclic (non‐ring) aliphatic molecule with exactly one terminal carboxylic acid group.
Furthermore, the sole carbon bonded to that acid group (the “alpha” carbon) must carry at least one hydroxyl (-OH)
substituent. To avoid misclassifications, we also require that the remainder of the molecule (the fatty acid chain)
be sufficiently “aliphatic”–namely, that the connected carbon chain extending from the alpha carbon is of a
minimum length.
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    Our criteria are:
      1) Molecule is acyclic (contains no rings).
      2) Molecule contains a terminal carboxylic acid group.
         - We match a carboxyl group via SMARTS: "C(=O)[O;H1]".
         - We further require that there is exactly one unique acid carbon (the carbon in the C=O).
      3) The single carbon attached to the acid carbon (the alpha carbon) must have a hydroxyl group.
      4) The alpha carbon must be part of an extended aliphatic chain (we check that the connected chain 
         of carbon atoms (ignoring the acid branch) is at least a minimal length).
        
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a 2-hydroxy fatty acid, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Reject molecules with rings – most fatty acids are acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a typical fatty acid"
    
    # 2. Look for terminal carboxylic acid group.
    # SMARTS for carboxylic acid: acidic carbon (C=O) attached to an -OH.
    carboxyl_smarts = "C(=O)[O;H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    # Get unique acid carbon indices from these matches (atom 0 in our SMARTS is the acid carbon)
    acid_carbons = set(match[0] for match in carboxyl_matches)
    if len(acid_carbons) != 1:
        return False, "Molecule must have exactly one terminal carboxylic acid group"
    
    acid_carbon_idx = list(acid_carbons)[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # 3. Identify the alpha carbon.
    # For a terminal carboxylic acid, the acid carbon should be attached to exactly one carbon (ignoring the oxygens).
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in list(acid_carbons)]
    if len(carbon_neighbors) != 1:
        return False, "Acid carbon is not terminal (attached to 0 or >1 carbons)"
    
    alpha_carbon = carbon_neighbors[0]
    
    # Check that alpha carbon has at least one hydroxyl substituent (an oxygen with at least one hydrogen)
    hydroxyl_found = False
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetIdx() == acid_carbon_idx:
            continue
        if nbr.GetAtomicNum() == 8:
            # total number of explicit hydrogens attached may be reported in GetTotalNumHs method.
            if nbr.GetTotalNumHs() > 0:
                hydroxyl_found = True
                break
    if not hydroxyl_found:
        return False, "Alpha carbon does not bear a hydroxyl (-OH) substituent"
    
    # 4. Check that the fatty acid chain is sufficiently long and aliphatic.
    # We traverse carbon atoms connected to the alpha carbon (excluding the acid branch) 
    # and count the number in the connected subgraph of carbons.
    def dfs(atom, visited):
        count = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                visited.add(nbr.GetIdx())
                count = max(count, 1 + dfs(nbr, visited))
        return count

    # Do not traverse back to the acid carbon.
    visited = {acid_carbon_idx, alpha_carbon.GetIdx()}
    chain_length = 1 + dfs(alpha_carbon, visited)  # count alpha carbon plus its extendable carbons
    # Heuristically require a chain length of at least 4 (this cutoff helps remove small diacids and polyfunctional acids)
    if chain_length < 4:
        return False, f"Carbon chain too short ({chain_length} carbons) for a fatty acid"
    
    return True, ("Found a terminal carboxylic acid with an alpha carbon bearing a hydroxyl group and a sufficiently long "
                  f"aliphatic chain ({chain_length} carbons), classifying the molecule as a 2-hydroxy fatty acid")

# Example usage:
if __name__ == "__main__":
    # Test with one of the positive examples: (2R)-2-hydroxytetradecanoic acid
    test_smiles = "CCCCCCCCCCCC[C@@H](O)C(O)=O"
    result, reason = is_2_hydroxy_fatty_acid(test_smiles)
    print(result, reason)