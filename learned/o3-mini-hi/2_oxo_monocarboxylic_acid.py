"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
A valid molecule must have exactly one carboxylic acid group (either as -C(=O)[OH] or -C(=O)[O-])
and the carbon atom directly attached to that carboxyl carbon (the α‐carbon) must itself have a double‐bonded oxygen (ketone functionality).
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if the given molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        (bool, str): A tuple where the boolean indicates if it is classified as a 2-oxo monocarboxylic acid 
                     and the string provides a detailed reasoning. If the SMILES is invalid returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, we want exactly one carboxylic acid group.
    # We use a SMARTS that covers both protonated and deprotonated forms.
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    
    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, expected exactly one (monocarboxylic acid)"
    
    # In our SMARTS match the first atom (index 0) is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # In a typical monocarboxylic acid the carboxyl carbon is terminal and hence should have exactly one carbon neighbor.
    alpha_candidates = []
    for neighbor in acid_carbon.GetNeighbors():
        # Exclude the oxygen(s) from the acid group; only count carbons.
        if neighbor.GetAtomicNum() == 6:
            alpha_candidates.append(neighbor)
    if len(alpha_candidates) == 0:
        return False, "No alpha carbon found (likely a formic acid derivative) so no possibility for a 2-oxo substituent"
    if len(alpha_candidates) > 1:
        return False, f"Found {len(alpha_candidates)} alpha carbon candidates – ambiguous for a 2-oxo substituent"
    
    alpha = alpha_candidates[0]
    # Now, require that the α‐carbon (directly attached to the carboxyl carbon) is also bonded (via a double bond) to an oxygen.
    found_2oxo = False
    for bond in alpha.GetBonds():
        # Look for a double bond.
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            other = bond.GetOtherAtom(alpha)
            # We expect the carbonyl oxygen to be attached (and ignore the case in which the candidate might be the acid carbon)
            if other.GetAtomicNum() == 8:
                found_2oxo = True
                break
    if not found_2oxo:
        return False, "Alpha carbon does not have a 2-oxo (carbonyl) substituent"
    
    return True, "Contains a single carboxylic acid group with an α‐carbon bearing a 2‑oxo (carbonyl) substituent"

# Example usage:
# Uncomment the following lines to test with an example SMILES.
# test_smiles = "CCCCCCC(=O)C(O)=O"  # 2-oxononanoic acid
# result, reason = is_2_oxo_monocarboxylic_acid(test_smiles)
# print(result, reason)