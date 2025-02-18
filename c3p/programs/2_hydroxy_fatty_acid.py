"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-hydroxy fatty acid
Definition: A fatty acid (an aliphatic molecule with a carboxylic acid group)
with a hydroxyl group on the alpha-carbon (i.e. the only carbon attached to the acid group).
The acid group must be terminal (attached to only one carbon) so that the “alpha” position is unambiguously defined.
Molecules with rings (e.g. sugars or macrocycles) or with the acid group embedded (acid carbon attached to >1 carbon) 
are not classified as fatty acids.
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    It checks that the molecule is acyclic, has a single (terminal) carboxylic acid group,
    and that the lone (alpha) carbon attached to that acid group carries a hydroxyl (-OH) substituent.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a 2-hydroxy fatty acid, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with rings (many fatty acids are acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a typical fatty acid"

    # Look for terminal carboxylic acid group.
    # SMARTS for carboxylic acid: a carbonyl carbon attached to an -OH.
    carboxyl_smarts = "C(=O)[O;H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Evaluate each matched carboxyl group.
    for match in carboxyl_matches:
        # In our SMARTS the first atom (index 0) is the acid carbon.
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        
        # For a terminal carboxylic acid, the acid carbon should be bound to exactly one carbon
        # that is not part of the carboxyl group (i.e. not one of the oxygens).
        carbon_neighbors = [
            neighbor for neighbor in acid_carbon.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in match
        ]
        if len(carbon_neighbors) != 1:
            # If the acid carbon is attached to 0 or >1 carbons (e.g., embedded in a ring or a polyacid),
            # then skip this match.
            continue
        
        # The only carbon neighbor is our candidate alpha carbon.
        alpha_carbon = carbon_neighbors[0]

        # Check the alpha carbon for a hydroxyl substituent.
        # We ignore the bond going back to the acid carbon.
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetIdx() == acid_carbon_idx:
                continue
            # We look for an oxygen that carries at least one hydrogen.
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                return True, ("Found terminal carboxylic acid with an alpha carbon bearing a hydroxyl group, "
                              "classifying the molecule as a 2-hydroxy fatty acid")
    # If none of the carboxylic acid groups satisfy the criteria…
    return False, "No terminal carboxylic acid with an alpha-hydroxy substituent found"

# Example usage:
if __name__ == "__main__":
    # (2R)-2-hydroxytetradecanoic acid (example of a true positive)
    test_smiles = "CCCCCCCCCCCC[C@@H](O)C(O)=O"
    result, reason = is_2_hydroxy_fatty_acid(test_smiles)
    print(result, reason)