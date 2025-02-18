"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone (acyloin)
Definition:
    A secondary α‐hydroxy ketone is an α‐hydroxy ketone in which the carbonyl group and the 
    hydroxy group are linked by a carbon bearing one hydrogen and one organyl group.
    This motif corresponds to R–CH(OH)–C(=O)–R′.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin) based on its SMILES string.
    The structure must contain the motif: R–CH(OH)–C(=O)–R′, where the CH bearing -OH is secondary.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains the secondary alpha-hydroxy ketone motif, False otherwise.
        str : Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the acyloin (secondary alpha-hydroxy ketone) motif.
    # Explanation of the SMARTS:
    #   [C;H1]: a carbon with exactly one hydrogen (making it secondary)
    #   ([OX2H]): bearing an -OH group (O atom with two connections and one hydrogen)
    #   C(=O): next is a carbonyl group
    #   [#6]: the carbonyl carbon is attached to an organic substituent (a carbon atom)
    acyloin_smarts = "[C;H1]([OX2H])C(=O)[#6]"
    acyloin_pattern = Chem.MolFromSmarts(acyloin_smarts)
    
    # Check whether the molecule has at least one match to the acyloin SMARTS pattern.
    if mol.HasSubstructMatch(acyloin_pattern):
        return True, "Molecule contains the secondary alpha-hydroxy ketone (acyloin) motif"
    else:
        return False, "Does not contain the required secondary alpha-hydroxy ketone motif"

# Optional: Test examples can be added here if run as a script. 
# For instance:
# if __name__ == '__main__':
#     test_smiles = "O[C@H](C(=O)c1ccccc1)c1ccccc1"  # (S)-benzoin example
#     result, reason = is_secondary_alpha_hydroxy_ketone(test_smiles)
#     print(result, reason)