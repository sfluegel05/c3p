"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
A valid molecule must have exactly one carboxylic acid group 
(–C(=O)[OH] or –C(=O)[O-]) and a single carbon neighbor (the α‐carbon)
that bears a carbonyl (C=O) substituent.
This improved version relaxes overly strict requirements on hybridization and connectivity.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    The function first verifies that there is exactly one carboxylic acid group.
    Then it finds the carbon (α‐carbon) directly bonded to the carboxyl carbon and checks 
    if this α‐carbon has at least one double-bond to an oxygen (the 2‑oxo substituent).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): Tuple with a boolean classification and a reasoning string.
                     Returns (False, "Invalid SMILES string") if the molecule cannot be parsed.
    """
    
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a carboxylic acid group (covers both -COOH and -COO- forms)
    acid_smarts = "[CX3](=O)[OX2H1,O-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    
    # Must have exactly one carboxylic acid group.
    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups; expected exactly one (monocarboxylic acid)"
    
    # In the SMARTS match, the first atom (index 0) is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Get all carbon neighbors of the acid carbon (potential α-carbons)
    alpha_candidates = []
    for neighbor in acid_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Only carbons
            alpha_candidates.append(neighbor)
    
    if not alpha_candidates:
        return False, "No alpha carbon found (likely a formic acid derivative) so no possibility for a 2-oxo substituent"
    
    valid_alpha = []
    # Check each α‐candidate: it must directly have a double bond to oxygen.
    for alpha in alpha_candidates:
        has_double_oxygen = False
        for bond in alpha.GetBonds():
            # Check for a double bond
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(alpha)
                if other.GetAtomicNum() == 8:  # Oxygen
                    has_double_oxygen = True
                    break
        if has_double_oxygen:
            valid_alpha.append(alpha)
    
    if len(valid_alpha) == 0:
        return False, "No α‐carbon with a C=O substituent found; molecule lacks a clear 2‑oxo group"
    if len(valid_alpha) > 1:
        return False, f"Found {len(valid_alpha)} potential α‐carbons with a C=O substituent — ambiguous for classification"
    
    return True, "Contains a single carboxylic acid group with an α‐carbon bearing a 2‑oxo (C=O) substituent"

# Example usage:
# Uncomment the lines below to test with one of the provided SMILES.
# test_smiles = "CCCCCCCC(=O)C(O)=O"  # 2-oxononanoic acid
# result, reason = is_2_oxo_monocarboxylic_acid(test_smiles)
# print(result, reason)