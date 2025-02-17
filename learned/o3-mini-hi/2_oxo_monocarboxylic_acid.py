"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
A valid molecule must have exactly one carboxylic acid group (either as -C(=O)[OH] or -C(=O)[O-])
and the carbon directly bonded to that acid carbon (the α‐carbon) must have a double bond to oxygen.
This program now includes additional checks on the alpha carbon’s connectivity and hybridization 
to reduce false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    The method first verifies that there is exactly one carboxylic acid group. Then it collects
    the carbon neighbors (alpha carbons) of the carboxyl carbon. For each candidate alpha carbon,
    it requires that:
      - The candidate is sp2-hybridized (as expected for a carbonyl carbon),
      - It has exactly three heavy-atom connections (one from the acid group, one to a substituent,
        and one to a double-bonded oxygen), and
      - It has exactly one double bond to an oxygen that is not part of the acid group.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): Tuple with the boolean classification and a reasoning string.
                     Returns (False, "Invalid SMILES string") for an invalid molecule.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find carboxylic acid group. This SMARTS covers protonated and deprotonated forms.
    acid_smarts = "[CX3](=O)[OX2H1,O-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Must have exactly one carboxylic acid group.
    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, expected exactly one (monocarboxylic acid)"
    
    # In the SMARTS match, the first atom (index 0) is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Identify carbon neighbors of the acid carbon (potential α-carbons).
    alpha_candidates = []
    for neighbor in acid_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # consider only carbons
            alpha_candidates.append(neighbor)
            
    if len(alpha_candidates) == 0:
        return False, "No alpha carbon found (likely a formic acid derivative) so no possibility for a 2-oxo substituent"
    
    valid_alpha = []
    # Check each candidate for a 2-oxo substituent.
    for alpha in alpha_candidates:
        # For clarity, get list of bonds from the alpha carbon.
        bonds = alpha.GetBonds()
        
        # Check that the candidate is sp2 hybridized (typical for a carbonyl carbon).
        if alpha.GetHybridization() != rdchem.HybridizationType.SP2:
            continue
        
        # Count heavy atom neighbors of the alpha carbon.
        # (A genuine α-keto carbon should have exactly three: one from acid carbon,
        # one from its C=O and one additional substituent.)
        heavy_neighbors = [nbr for nbr in alpha.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 3:
            continue
        
        # Ensure one of the neighbors is the acid carbon.
        if acid_carbon not in heavy_neighbors:
            continue
        
        # Count double bonds from alpha to oxygen (that are not the acid group bond).
        dbl_oxygen_count = 0
        for bond in bonds:
            # Look for double bonds
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(alpha)
                # Skip if this oxygen is the one attached to the acid group (if accidentally shared)
                if other.GetAtomicNum() == 8 and other.GetIdx() != acid_carbon_idx:
                    dbl_oxygen_count += 1
        if dbl_oxygen_count == 1:
            valid_alpha.append(alpha)
    
    # We require exactly one unambiguous alpha candidate with the requisite carbonyl.
    if len(valid_alpha) == 0:
        return False, "Alpha carbon does not have a clear 2-oxo (carbonyl) substituent based on connectivity/hybridization"
    if len(valid_alpha) > 1:
        return False, f"Found {len(valid_alpha)} alpha carbon candidates with a 2-oxo substituent – ambiguous for classification"
    
    return True, "Contains a single carboxylic acid group with an α‐carbon bearing a 2‑oxo (carbonyl) substituent"

# Example usage:
# Uncomment the following lines to test an example SMILES:
# test_smiles = "CCCCCCCC(=O)C(O)=O"  # 2-oxononanoic acid
# result, reason = is_2_oxo_monocarboxylic_acid(test_smiles)
# print(result, reason)