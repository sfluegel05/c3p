"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
at the 1-position of the glycerol fragment.
Our improved heuristic checks two conditions:
  1. Presence of a phosphoethanolamine headgroup. We look for a substructure
     typical for glycerophosphoethanolamines: a phosphate connected to an ethanolamine.
     Here we use the SMARTS "COP(=O)(O)OCCN" (ignoring charges and stereochemistry).
  2. Presence of a typical glycerol fragment that carries an acyl (ester) substituent at the sn-1 position.
     Many of the known examples contain a fragment like "OCC(O)COC(=O)" (i.e. a glycerol where the sn-1 OH is esterified).
     
These two combined requirements help reduce false positives.
"""

from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    We require that the molecule satisfies:
      (a) It contains a phosphoethanolamine headgroup.
          We search for a substructure based on "COP(=O)(O)OCCN".
      (b) It contains a glycerol-like backbone with an esterified sn-1 hydroxyl.
          Typical examples show the motif "OCC(O)COC(=O)" where the terminal "OC(=O)" is the 1-O-acyl substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to be a 1-O-acylglycerophosphoethanolamine, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Condition 1: Check for the phosphoethanolamine headgroup.
    # This pattern looks for a phosphate (P) with one =O, one -OH, and one -OCCN substituent.
    headgroup_smarts = "COP(=O)(O)OCCN"
    headgroup_pattern = Chem.MolFromSmarts(headgroup_smarts)
    if headgroup_pattern is None:
        return False, "Error in headgroup SMARTS pattern"
    if not mol.HasSubstructMatch(headgroup_pattern):
        return False, "Phosphoethanolamine headgroup not found"
    
    # Condition 2: Check for the glycerol backbone with a 1-O-acyl substitution.
    # Many true examples include a fragment like "OCC(O)COC(=O)",
    # which represents a glycerol with the sn-1 hydroxyl esterified (OC(=O)...).
    glycerol_acyl_smarts = "OCC(O)COC(=O)"
    glycerol_acyl_pattern = Chem.MolFromSmarts(glycerol_acyl_smarts)
    if glycerol_acyl_pattern is None:
        return False, "Error in glycerol-acyl SMARTS pattern"
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "Glycerol backbone with 1-O-acyl substitution not clearly identified"
    
    # Optional extra check: You could add tests such as verifying the acyl chain is sufficiently long.
    # For example, one might check that the carbon chain attached to the C(=O) is at least 6 carbons.
    # Here we simply note that the presence of the glycerol-acyl pattern is the key criterion.
    
    return True, "Molecule contains a phosphoethanolamine headgroup with 1-O-acyl substitution on glycerol"

# Example usage:
# result, reason = is_1_O_acylglycerophosphoethanolamine("CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN")
# print(result, reason)