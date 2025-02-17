"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone (acyloin)

Definition:
    A secondary α‐hydroxy ketone (acyloin) is a compound containing the motif 
        R–CH(OH)–C(=O)–R′
    in which the CH bearing the –OH is secondary, meaning it is bonded to three heavy atoms:
    one from the hydroxyl group, one from the carbonyl group (C=O), and one organyl (organic) group.
    
Our previous attempt used a simple SMARTS,
    "[C;H1]([OX2H])C(=O)[#6]"
but this proved too unspecific. The updated SMARTS explicitly requires a carbon 
substituent on the left and that the acyloin center has exactly three heavy-atom neighbors.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone (acyloin) based on its SMILES string.
    The structure must contain the motif: R–CH(OH)–C(=O)–R′ where the CH bearing -OH is secondary.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains the secondary alpha-hydroxy ketone motif, False otherwise.
        str : Reason for classification.
    """
    # Parse the SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern:
    # 1. "[#6]-" specifies that the acyloin carbon has a carbon substituent (an organic group) on one side.
    # 2. "[C;X3;H1]([OX2H])" enforces that the acyloin carbon is tetrahedral with exactly three heavy atoms (making it secondary) 
    #    and that it bears an -OH group.
    # 3. "C(=O)[#6]" requires a carbonyl group next, where the carbonyl carbon is attached to another carbon.
    acyloin_smarts = "[#6]-[C;X3;H1]([OX2H])C(=O)[#6]"
    acyloin_pattern = Chem.MolFromSmarts(acyloin_smarts)
    
    if acyloin_pattern is None:
        return False, "Error processing SMARTS pattern"
    
    # Check for at least one match to the acyloin motif.
    matches = mol.GetSubstructMatches(acyloin_pattern)
    if matches:
        # Optionally, one can perform further checks per match.
        return True, "Molecule contains the secondary alpha-hydroxy ketone (acyloin) motif"
    else:
        return False, "Does not contain the required secondary alpha-hydroxy ketone motif"

# Optional: execute tests when run as a script.
if __name__ == '__main__':
    # Example test case: (S)-benzoin should be classified as an acyloin.
    test_smiles = "O[C@H](C(=O)c1ccccc1)c1ccccc1"
    result, reason = is_secondary_alpha_hydroxy_ketone(test_smiles)
    print(result, reason)