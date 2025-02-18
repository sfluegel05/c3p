"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 1,2-diacyl-sn-glycero-3-phosphocholine
Defined as: The conjugate base of a 1,2-diacyl-sn-glycero-3-phosphocholine compound 
formed by deprotonation of the phosphate OH group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.
    Expected features:
      1. A phosphocholine head group, recognizable by a phosphate group with a negative oxygen and
         an attached choline fragment (OCC[N+](C)(C)C).
      2. Exactly two acyl ester groups (OC(=O)) that represent the fatty acid chains linked to the glycerol backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule matches the class, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the phosphocholine head group.
    # This pattern looks for a phosphorus atom with:
    #   - a double bond O (P=O)
    #   - an oxygen with a negative charge ([O-])
    #   - an oxygen connected to an ethyl chain ending in a quaternary ammonium (OCC[N+](C)(C)C)
    pc_pattern = Chem.MolFromSmarts("P(=O)([O-])(OCC[N+](C)(C)C)")
    if not mol.HasSubstructMatch(pc_pattern):
        return False, "Phosphocholine head group not found"
    
    # Define a SMARTS pattern for an acyl ester group (fatty acid ester, part of the glycerol backbone).
    # This is the typical fragment "OC(=O)" linking a fatty acid chain.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # We expect exactly two acyl ester groups (for sn-1 and sn-2 acylation)
    if len(ester_matches) != 2:
        return False, f"Expected 2 acyl ester groups, found {len(ester_matches)}"
    
    # Optionally, one can perform additional checks such as ensuring that the molecule contains a glycerol-like backbone 
    # and verifying overall molecular properties, but for our present classification this suffices.
    
    return True, ("Molecule contains a phosphocholine head group with a deprotonated phosphate "
                  "and exactly two acyl ester groups, consistent with a 1,2-diacyl-sn-glycero-3-phosphocholine.")

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCCCCCCC)COC(=O)CCCCCCCCC)([O-])=O"
    result, reason = is_1_2_diacyl_sn_glycero_3_phosphocholine(test_smiles)
    print(result, reason)