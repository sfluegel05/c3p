"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
A valid molecule must have exactly one carboxylic acid group (either as -C(=O)[OH] or -C(=O)[O-])
and the carbon directly bonded to that acid carbon (the α‐carbon) must also be bonded (via a double bond)
to an oxygen (i.e. have a ketone functionality). This program implements additional checks to avoid ambiguity.
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    The method first verifies that there is exactly one carboxylic acid group. Then it collects
    the carbon neighbors (alpha carbons) of the carboxyl carbon. Finally, it checks that exactly one
    of these alpha carbons has at least one double bond to an oxygen (other than the acid group oxygen).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): Tuple with the boolean classification and a reasoning string.
                     Returns (False, "Invalid SMILES string") for an invalid molecule.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find carboxylic acid group. This SMARTS covers the protonated and deprotonated forms.
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
    
    # Identify all carbon neighbors of the carboxyl carbon.
    alpha_candidates = []
    for neighbor in acid_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Only consider carbons
            alpha_candidates.append(neighbor)
    
    if len(alpha_candidates) == 0:
        # In formic acid derivatives there is no alpha carbon.
        return False, "No alpha carbon found (likely a formic acid derivative) so no possibility for a 2-oxo substituent"
    
    # Now search among the alpha candidates for one that is bonded to an oxygen with a double bond.
    valid_alpha = []
    for alpha in alpha_candidates:
        for bond in alpha.GetBonds():
            # Ensure the bond is a double bond
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(alpha)
                # We ignore if the oxygen is the one already in the acid group.
                if other.GetAtomicNum() == 8 and other.GetIdx() != acid_carbon_idx:
                    valid_alpha.append(alpha)
                    break  # Stop after the first valid double bond is found for this candidate
    
    # We require exactly one unambiguous alpha candidate with the carbonyl substituent.
    if len(valid_alpha) == 0:
        return False, "Alpha carbon does not have a 2-oxo (carbonyl) substituent"
    if len(valid_alpha) > 1:
        return False, f"Found {len(valid_alpha)} alpha carbon candidates with a 2-oxo substituent – ambiguous for classification"
    
    # If we get here, we have exactly one carboxylic acid and one corresponding alpha carbon
    # with a double-bonded oxygen that is not part of the acid group.
    return True, "Contains a single carboxylic acid group with an α‐carbon bearing a 2‑oxo (carbonyl) substituent"

# Example usage:
# Uncomment the following lines to test an example SMILES.
# test_smiles = "CCCCCCC(=O)C(O)=O"  # 2-oxononanoic acid
# result, reason = is_2_oxo_monocarboxylic_acid(test_smiles)
# print(result, reason)