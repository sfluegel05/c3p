"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: CHEBI:1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
at the 1-position of the glycerol fragment.
Heuristic:
  We search for an integrated substructure pattern that encodes the expected backbone:
    - An acyl ester fragment: C(=O)O
    - Attached to a glycerol backbone: CC(O)CO (which stands in for CH2–CH(OH)–CH2)
    - Followed by a phosphate group bound to ethanolamine: [P](=O)(O)OCCN
  The overall SMARTS used is:
       "C(=O)OCC(O)CO[P](=O)(O)OCCN"
  This pattern is intended to capture that the acyl group is directly attached (via an ester bond)
  to the sn-1 oxygen of a glycerol fragment whose sn-3 oxygen is esterified to a phosphate that carries an ethanolamine.
  
  Note: Because many false positives arose from imprecise headgroup definitions (e.g.
        phosphoserine, phosphocholine), by requiring “OCCN” exactly in the pattern and the complete glycerol backbone,
        we reduce misclassifications.
"""

from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    The method uses an improved integrated SMARTS pattern:
      "C(=O)OCC(O)CO[P](=O)(O)OCCN"
    which encodes:
      - an acyl carbonyl (C(=O)) attached via an ester oxygen to a glycerol fragment (CC(O)CO)
      - the glycerol is further connected at sn-3 via an O to a phosphate [P](=O)(O)OCCN headgroup,
        where the headgroup is defined as an ethanolamine (OCCN).
        
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 1-O-acylglycerophosphoethanolamine, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the integrated SMARTS pattern.
    # This pattern looks for:
    #   Acyl group:                C(=O)O
    #   Glycerol backbone:          CC(O)CO    (representing CH2–CH(OH)–CH2; chirality is not enforced)
    #   Phosphoethanolamine head:   [P](=O)(O)OCCN 
    integrated_smarts = "C(=O)OCC(O)CO[P](=O)(O)OCCN"
    
    integrated_pattern = Chem.MolFromSmarts(integrated_smarts)
    if integrated_pattern is None:
        return False, "Error in constructing integrated SMARTS pattern"
    
    # Check for a substructure match and return appropriate result.
    if mol.HasSubstructMatch(integrated_pattern):
        return True, ("Molecule contains a phosphoethanolamine headgroup with a 1-O-acyl substitution "
                      "on glycerol; the acyl (carbonyl) group is directly attached via an ester at the sn-1 position.")
    else:
        return False, ("No matching integrated structure found. The molecule may lack the required 1-O-acyl substitution "
                       "on the glycerol fragment with a phosphoethanolamine headgroup.")

# Example usage:
# Uncomment the lines below to test one example.
# test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
# result, explanation = is_1_O_acylglycerophosphoethanolamine(test_smiles)
# print(result, explanation)