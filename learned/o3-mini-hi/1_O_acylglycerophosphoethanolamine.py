"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
at the 1-position of the glycerol fragment.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    For classification we check:
      (1) The presence of a phosphoethanolamine headgroup. 
          We use the SMARTS pattern "COP(=O)(O)OCCN" which looks for a phosphate (P with =O and an -OCCN substituent).
      (2) The presence of an ester linkage at a primary alcohol group.
          Here we use "[CH2]OC(=O)" to suggest an acyl group bound via an ester to a CH2 group.
      (3) A glycerol-like backbone (heuristic) because the glycerol fragment typically would convey connectivity.
    
    These patterns are heuristic and are intended to capture examples such as:
      - 1-arachidonoyl-sn-glycero-3-phosphoethanolamine
      - 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
      - 1-oleoyl-sn-glycero-3-phosphoethanolamine, etc.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a 1-O-acylglycerophosphoethanolamine, False otherwise.
        str: Reason explaining the classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the phosphoethanolamine headgroup.
    # This pattern looks for an oxygen directly attached to a phosphorus which has =O, one free -OH (or equivalent)
    # and an -OCCN moiety.
    headgroup_smarts = "COP(=O)(O)OCCN"
    headgroup_pattern = Chem.MolFromSmarts(headgroup_smarts)
    if not mol.HasSubstructMatch(headgroup_pattern):
        return False, "Phosphoethanolamine headgroup not found"
    
    # Check for the ester linkage at the 1-position.
    # This pattern "[CH2]OC(=O)" is meant to mark an ester bond where an acyl (C(=O)R) is linked via oxygen 
    # to a CH2 unit (expected from the sn-1 position of glycerol).
    acyl_smarts = "[CH2]OC(=O)"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "1-O-acyl (ester linkage) substitution on glycerol not found"
    
    # (Optional) Check for a glycerol-like backbone.
    # For a 1-O-acylglycerophosphoethanolamine, the glycerol fragment would ideally show a 3-carbon motif.
    # Since one of the primary –OH groups is now esterified we use a simplified pattern.
    glycerol_smarts = "O[C;H2]C(O)C"  # This loosely captures a fragment with a free –OH in the central carbon.
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not clearly identified"
    
    return True, "Molecule contains a phosphoethanolamine headgroup with 1-O-acyl substitution on glycerol."

# For example, to test the function you might call:
# result, reason = is_1_O_acylglycerophosphoethanolamine("CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN")
# print(result, reason)