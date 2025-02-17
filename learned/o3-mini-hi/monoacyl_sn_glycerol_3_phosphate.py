"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate
Definition:
  An sn-glycero-3-phosphate compound having a single unspecified acyl group at either position 1 or position 2.
  
Approach:
  1. Parse the SMILES.
  2. Check that there is exactly one phosphorus atom and that a phosphate substructure is present.
  3. Identify an acyl ester moiety ([CX3](=O)O[CX4]) that is not a phosphate ester (i.e. the oxygen is not attached to P).
     There must be exactly one such acyl ester.
  4. Look for a simplified glycerol-phosphate backbone pattern.
  
If any of these tests fail, the function returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if the molecule given by the SMILES string is a monoacyl-sn-glycerol 3-phosphate.
    
    A monoacyl-sn-glycerol 3-phosphate is defined as an sn-glycero-3-phosphate compound having 
    a single acyl group (an esterified fatty acid) at either position 1 or position 2.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES; return error if parsing fails.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for exactly one phosphorus (the phosphate group) ---
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Molecule must contain exactly one phosphorus atom (found {len(phosphorus_atoms)})"
    
    # --- Check for a phosphate group ---
    # This SMARTS checks for a phosphorus with one double-bonded oxygen and three oxygens.
    phosphate_smarts = "P(=O)(O)(O)"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # --- Identify the acyl ester group ---
    # We look for an ester linkage: a carbonyl (C(=O)) bound to O which is bound to an alkyl carbon.
    # The SMARTS [CX3](=O)O[#6] should match a typical acyl ester fragment.
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)O[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    acyl_count = 0
    for match in acyl_matches:
        # match is a tuple of atom indices: (carbonyl carbon, ester oxygen, alkyl carbon)
        o_atom = mol.GetAtomWithIdx(match[1])
        # Exclude ester groups where the oxygen is attached to a phosphorus
        if any(neighbor.GetAtomicNum() == 15 for neighbor in o_atom.GetNeighbors()):
            continue
        acyl_count += 1
    if acyl_count != 1:
        return False, f"Found {acyl_count} acyl ester group(s); expected exactly one for a monoacyl-sn-glycerol 3-phosphate"
    
    # --- Check for a glycerol 3-phosphate backbone ---
    # We use a simplified pattern for a glycerol fragment bearing a phosphate:
    #   A primary (CH2) carbon, a secondary (CH) carbon, and a third (CH2) carbon bound to a phosphate
    # This pattern is an approximation.
    glycerol_phosphate_smarts = "[CX4;H2][CX4;H][CX4;H2](OP(=O)(O)O)"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol 3-phosphate backbone not found"
    
    return True, "Molecule has a glycerol 3-phosphate backbone with exactly one acyl ester group"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one of the provided SMILES strings:
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O"
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, reason)