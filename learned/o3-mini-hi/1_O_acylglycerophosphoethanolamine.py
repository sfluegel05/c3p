"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
at the 1-position of the glycerol fragment.
Heuristic:
  Instead of first matching the GPE headgroup and then verifying the acyl substitution,
  a combined SMARTS query is employed:
      C(=O)OCC(O)COP(=O)(O)OCCN
  This pattern encodes:
    - An acyl carbonyl (C(=O))
    - Followed by an ester oxygen linked to a glycerol fragment (OCC(O)CO) 
      where the glycerol is represented as CH2–CH(OH)–CH2.
    - Followed by the phosphate and ethanolamine headgroup (P(=O)(O)OCCN).
  A match indicates that the acyl group is directly attached (via an ester bond)
  at the sn-1 position of a glycerophosphoethanolamine.
"""

from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as a 1-O-acylglycerophosphoethanolamine, else False.
      str: Explanation of the decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Construct an integrated SMARTS pattern which encodes:
    #  - an acyl carbonyl group (C(=O))
    #  - an ester oxygen immediately attached to a glycerol fragment (OCC(O)CO)
    #  - the glycerol is then linked via the sn3 position to a phosphate that carries an ethanolamine (OCCN)
    # If the molecule contains this unit, we assume that it is a correct 1-O-acylglycerophosphoethanolamine.
    integrated_smarts = "C(=O)OCC(O)COP(=O)(O)OCCN"
    integrated_pattern = Chem.MolFromSmarts(integrated_smarts)
    if integrated_pattern is None:
        return False, "Error in constructing integrated SMARTS pattern"

    # Check for a substructure match
    if mol.HasSubstructMatch(integrated_pattern):
        return True, ("Molecule contains a phosphoethanolamine headgroup with a 1-O-acyl substitution "
                      "on glycerol; the acyl group (carbonyl) is directly connected to the sn-1 oxygen.")
    else:
        return False, ("No matching integrated structure found. The molecule may lack the required "
                       "1-O-acyl substitution on the glycerol fragment with the expected phosphoethanolamine headgroup.")

# Example usage:
# Uncomment the following lines to test one of the provided SMILES.
# smiles_example = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"  # 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine
# result, reason = is_1_O_acylglycerophosphoethanolamine(smiles_example)
# print(result, reason)