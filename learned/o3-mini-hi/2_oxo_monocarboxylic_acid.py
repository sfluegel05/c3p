"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-oxo monocarboxylic acid
Definition: Any monocarboxylic acid having a 2-oxo substituent.
The molecule must contain exactly one carboxylic acid group (-C(=O)O or its deprotonated form)
and the carbon atom adjacent (α‐carbon) to the carboxyl carbon must be attached to a carbonyl (C=O) bond.
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if the given molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        (bool, str): Tuple with a boolean indicating if it is a 2-oxo monocarboxylic acid
                     and a string with the reasoning.
                     If the input SMILES is invalid, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for carboxylic acid.
    # One pattern for the protonated acid and one for the deprotonated form.
    acid_smarts1 = Chem.MolFromSmarts("C(=O)[OX2H]")
    acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")
    
    # Find matches for each pattern.
    acid_matches1 = mol.GetSubstructMatches(acid_smarts1)
    acid_matches2 = mol.GetSubstructMatches(acid_smarts2)
    
    # Combine matches (each match is a tuple of atom indices where index 0 is the carboxyl carbon).
    acid_matches = acid_matches1 + acid_matches2

    if len(acid_matches) == 0:
        return False, "No carboxylic acid group found"
    if len(acid_matches) > 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, expected exactly one (monocarboxylic acid)"
    
    # Use the unique match.
    # In our SMARTS, the first atom (index 0) corresponds to the acid carbon (C=O part)
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Find the alpha carbon: it is a carbon neighbor of the acid carbon that is not the carbonyl oxygen.
    alpha_carbons = []
    for neighbor in acid_carbon.GetNeighbors():
        # Only consider carbon neighbors.
        if neighbor.GetAtomicNum() == 6:
            alpha_carbons.append(neighbor)
    if len(alpha_carbons) == 0:
        return False, "No alpha carbon found (likely formic acid), so no possibility for a 2-oxo substituent"
    
    # For each candidate alpha carbon, check if it is attached via a double bond to an oxygen (ketone functionality)
    found_2oxo = False
    for alpha in alpha_carbons:
        # Loop over bonds from the alpha carbon.
        for bond in alpha.GetBonds():
            # Check if the bond is a double bond.
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                # Get the atom on the other end.
                other = bond.GetOtherAtom(alpha)
                # Check if the other atom is oxygen.
                if other.GetAtomicNum() == 8:
                    # Found a carbonyl group attached to the alpha carbon.
                    found_2oxo = True
                    break
        if found_2oxo:
            break

    if not found_2oxo:
        return False, "Alpha carbon does not have a 2-oxo (carbonyl) substituent"

    return True, "Contains a single carboxylic acid group with an alpha carbon having a carbonyl (2-oxo) substituent"

# Example usage: Uncomment the below lines to test the function with one of the provided SMILES.
# test_smiles = "CCCCCCCC(=O)C(O)=O"  # 2-oxononanoic acid
# result, reason = is_2_oxo_monocarboxylic_acid(test_smiles)
# print(result, reason)