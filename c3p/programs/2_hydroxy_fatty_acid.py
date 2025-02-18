"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-hydroxy fatty acid
Definition: Any fatty acid with a hydroxy functional group in the alpha- or 2-position.
A fatty acid is essentially a molecule with a carboxylic acid group and an aliphatic chain.
This function checks for a carboxylic acid group and then verifies that the alpha carbon
(has a bond to the carboxyl carbon) is hydroxylated.
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    It looks for a carboxylic acid (C(=O)O) group and then inspects the carbon directly
    attached to the acid group (the alpha carbon) for a hydroxyl (–OH) group.
    
    Args:
        smiles (str): A SMILES string representing the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a 2-hydroxy fatty acid, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse SMILES and check validity
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group in the molecule.
    # SMARTS for carboxylic acid: a carbonyl carbon attached to an -OH.
    # Note: [O;H1] means oxygen with one hydrogen.
    carboxyl_smarts = "C(=O)[O;H1]"
    carboxyl_group = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_group)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    # For each carboxyl group found, check the alpha carbon (the neighbor of the acid carbon
    # that is not the oxygens of the acid group) for a hydroxyl (-OH) group.
    for match in carboxyl_matches:
        # In the SMARTS match, the first atom (index 0 in the match) is the acid carbon.
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        
        # Find neighbors of the carboxyl carbon that are not part of the carboxyl group oxygens.
        neighbor_indices = []
        for neighbor in acid_carbon.GetNeighbors():
            if neighbor.GetIdx() not in (match[1], match[2]):
                if neighbor.GetAtomicNum() == 6:  # check that neighbor is carbon
                    neighbor_indices.append(neighbor.GetIdx())
        
        if not neighbor_indices:
            # This carboxyl group might be terminally attached (no alpha carbon found)
            continue
        
        # For each candidate alpha carbon, check if it carries an -OH group.
        for alpha_idx in neighbor_indices:
            alpha_carbon = mol.GetAtomWithIdx(alpha_idx)
            # Check neighbors of the alpha carbon (excluding the carboxyl carbon)
            for atom in alpha_carbon.GetNeighbors():
                if atom.GetIdx() == acid_carbon_idx:
                    continue
                # The atom is a candidate for a hydroxyl oxygen
                if atom.GetAtomicNum() == 8:
                    # Verify that the oxygen has at least one hydrogen.
                    # GetTotalNumHs considers implicit hydrogens as well.
                    if atom.GetTotalNumHs() > 0:
                        return True, ("Found a carboxylic acid group with an alpha carbon bearing a hydroxyl "
                                      "group, classifying the molecule as a 2-hydroxy fatty acid")
    
    return False, "No alpha-hydroxy group found on any carboxylic acid group (2-position)"
    
# Example usage
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCC[C@@H](O)C(O)=O"  # (2R)-2-hydroxytetradecanoic acid
    result, reason = is_2_hydroxy_fatty_acid(test_smiles)
    print(result, reason)