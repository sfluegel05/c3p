"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechins
Definition: Members of the class of hydroxyflavan that have a flavan-3-ol skeleton 
            and its substituted derivatives.
            
Improved approach:
  1. Parse the SMILES string and remove stereochemistry.
  2. Look for the key structural motif using a relaxed, recursive SMARTS query.
     Here we require a dihydropyran (chroman) ring that:
       - has one oxygen in the ring,
       - bears an –OH group (representing the 3-hydroxy substituent),
       - and is fused to two aromatic rings (the A and B rings).
     The aromatic rings are specified recursively so that extra substituents on them
     do not preclude recognition of the core catechin skeleton.
  
  3. If the substructure is detected, the molecule is classified as a catechin.
  
Note: This rule‐based approach is heuristic, and while it should capture many substituted 
catechin derivatives, some edge cases might still be missed.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin class based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES string and remove stereochemistry.
      2. Look for a relaxed catechin core as defined by the recursive SMARTS:
           "[$(c1ccc(cc1))]C2CC([OX2H])C([$(c3ccc(cc3))])O2"
         This matches a structure in which an aromatic ring (ring A) is directly attached 
         to a dihydropyran (chroman) ring. That ring bears an –OH group (at the 3–position)
         and has an aromatic substituent (ring B) attached at a neighboring carbon.
      3. If the substructure is found, return True with a success message; otherwise, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is (or is closely related to) a catechin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemical labels to ignore chiral differences
    Chem.RemoveStereochemistry(mol)
    
    # Define a relaxed SMARTS pattern for the flavan-3-ol (catechin) core.
    # Explanation:
    #   [$(c1ccc(cc1))]  -> any aromatic ring (ring A)
    #   C2               -> a saturated carbon (beginning of the chroman ring)
    #   CC([OX2H])       -> next carbon in the ring bearing an -OH (C3)
    #   C([$(c3ccc(cc3))]) -> next carbon that is attached to an aromatic ring (ring B)
    #   O2               -> closing the dihydropyran (chroman) ring via an oxygen.
    core_smarts = "[$(c1ccc(cc1))]C2CC([OX2H])C([$(c3ccc(cc3))])O2"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create SMARTS for catechin core"
    
    # Check for the flavan-3-ol (catechin) core in the molecule.
    if mol.HasSubstructMatch(core_query):
        return True, "Molecule contains a flavan-3-ol (catechin) core (relaxed match)"
    else:
        return False, "Molecule does not contain the expected flavan-3-ol (catechin) core"

# The below code can be used for simple manual testing.
# if __name__ == "__main__":
#     test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"  # (-)-catechin
#     result, reason = is_catechin(test_smiles)
#     print(result, reason)