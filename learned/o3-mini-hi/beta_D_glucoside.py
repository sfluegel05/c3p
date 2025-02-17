"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: beta-D-glucoside
Definition: Any D-glucoside in which the anomeric centre has beta-configuration.
In our approach we attempt to detect the characteristic beta-D-glucopyranoside unit.
A typical unmodified beta-D-glucopyranoside fragment is illustrated by the SMARTS pattern:
    [!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
Here:
  • [!H0] stands for a non-hydrogen atom (from the aglycone).
  • The bond “-O” is the glycosidic oxygen.
  • [C@@H]1 is the anomeric carbon drawn with @@ stereochemistry (interpreted as beta).
  • The six-membered ring (positions 1–5 plus the ring oxygen) is built with typical substituents:
      – A CH2OH (CO) at one position, and –OH groups at the others.
Because sometimes the CH2OH group is modified (e.g. acetylated) we also allow a variant 
where “CO” is replaced by C(O*) (i.e. a carbon bearing an oxygen substituent, with wildcard).
"""

from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    We define a beta-D-glucoside as a molecule containing a pyranose (six-membered sugar) ring
    joined via an oxygen to an aglycone and having the anomeric carbon drawn with the beta (i.e.
    "@@") stereochemistry.
    
    We use two SMARTS substructure patterns:
      Pattern A (for the unmodified sugar):
          [!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O
      Pattern B (allowing modification at the CH2OH group):
          [!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a beta-D-glucoside, False otherwise.
        str : Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for beta-D-glucosides.
    # Pattern A: exactly matches common beta-D-glucopyranoside (no substitution on the CH2OH group)
    patternA_smarts = "[!H0]-O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    # Pattern B: a more flexible pattern allowing the CH2OH group to be substituted (e.g. acetylated)
    patternB_smarts = "[!H0]-O[C@@H]1O[C@H](C(O*))[C@@H](O)[C@H](O)[C@H]1O"
    
    patternA = Chem.MolFromSmarts(patternA_smarts)
    patternB = Chem.MolFromSmarts(patternB_smarts)
    
    if patternA is None or patternB is None:
        return False, "Error in SMARTS pattern definition"
    
    # Check if either of the patterns is found
    if mol.HasSubstructMatch(patternA) or mol.HasSubstructMatch(patternB):
        return True, "Molecule contains a beta-D-glucoside moiety with beta-anomeric configuration."
    else:
        return False, "No beta-D-glucoside substructure (with beta anomeric centre) found."

# Example test (uncomment to run)
# test_smiles = [
#     "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # beta-D-glucosyl-N-eicosanoylsphingosine
#     "[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C"  # oleuropein
# ]
# for s in test_smiles:
#     res, reason = is_beta_D_glucoside(s)
#     print(s, "=>", res, reason)