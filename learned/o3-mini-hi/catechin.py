"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechins
Definition: Members of the class of hydroxyflavan that have a flavan-3-ol skeleton 
            and its substituted derivatives.
            
Improved approach:
  1. Parse the SMILES and remove any stereochemical labels to simplify matching.
  2. Look for an explicit flavan-3-ol core using a SMARTS pattern that encodes
     a 2-phenyl-3,4-dihydro-2H-chromene core with an OH (or its substituted variant)
     at the 3-position. 
     
     The SMARTS pattern used is:
       "c1ccc(c(c1))C2CC(O)C(c3ccccc3)O2"
       
     This pattern requires an aromatic ring fused to a dihydropyran ring that bears 
     an –OH at one of its saturated carbons and that carries an aromatic substituent.
     
  3. If this core is detected the molecule is classified as a catechin.
     
Note:
   This rule-based approach is heuristic. The SMARTS pattern is strict so that only 
   molecules with an intact flavan-3-ol skeleton are returned; many derivatives may 
   have additional substituents, but as long as the core connectivity is maintained the 
   molecule is classified as catechin.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin class based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES string.
      2. Remove stereochemistry to simplify the pattern matching.
      3. Look for the flavan-3-ol core defined by the SMARTS:
            "c1ccc(c(c1))C2CC(O)C(c3ccccc3)O2"
         This represents a 2-phenyl-3,4-dihydro-2H-chromene core,
         with an -OH attached to the saturated C (at position 3) and a phenyl group (ring B)
         attached at the adjacent carbon.
      4. Return True with an appropriate message if the substructure is found, otherwise False.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains the flavan-3-ol core, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemistry to avoid mismatches due to chiral labels.
    Chem.RemoveStereochemistry(mol)
    
    # Define a SMARTS for the flavan-3-ol (catechin) core.
    # The pattern encodes:
    # - An aromatic ring (ring A) fused to a dihydropyran ring (ring C)
    # - The dihydropyran ring carries an -OH at one of its saturated carbons (the 3–position)
    # - An aromatic substituent (ring B) is attached to the adjacent saturated carbon.
    core_smarts = "c1ccc(c(c1))C2CC(O)C(c3ccccc3)O2"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create SMARTS pattern for catechin core"
    
    # Search for the catechin core substructure in the molecule.
    if not mol.HasSubstructMatch(core_query):
        return False, "Molecule does not contain the expected flavan-3-ol (catechin) core"
    
    return True, "Molecule contains a flavan-3-ol (catechin) core"

# Uncomment the following lines to run a simple test:
# if __name__ == "__main__":
#     test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"  # (-)-catechin
#     result, reason = is_catechin(test_smiles)
#     print(result, reason)