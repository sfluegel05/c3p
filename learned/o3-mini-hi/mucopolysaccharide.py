"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: Mucopolysaccharide (glycosaminoglycan)
Definition: “Any of the group of polysaccharides composed of alternating units from uronic acids and glycosamines,
and commonly partially esterified with sulfuric acid.”
Note:
  Because mucopolysaccharides are polymers made of sugar‐derived units, it is challenging to classify them
  solely from a SMILES string. This code uses heuristic SMARTS patterns to look for “uronic acid” and “glycosamine”
  rings (in a fully unconstrained “non‐stereospecific” manner) and reports True if a good number of each is found and roughly alternating.
  If the problem is too complex to solve reliably, the function may return (None, None).
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    The heuristic criteria are:
      1. The molecule must be parsed successfully.
      2. It should contain multiple sugar–like rings.
      3. At least two rings matching a uronic acid motif are found.
         (Here we use a simplified SMARTS pattern for a six‐membered ring with a carboxylic acid group.)
      4. At least two rings matching a glycosamine motif are found.
         (Here we use a simplified SMARTS pattern for a six–membered ring with an amino group.)
      5. The counts of these two types of rings should be roughly equal (to approximate an alternating pattern).
      6. If sulfate ester groups (–OS(=O)(=O)O) are present, that is noted.
    
    Note: This analysis is a rough heuristic and may mis‐classify some compounds.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as a mucopolysaccharide, False otherwise.
      str: Reason for the classification decision.
    """
    # First, attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define heuristic SMARTS for a uronic acid unit.
    # This is a simplified pattern for a pyranose ring bearing a carboxylic acid.
    # For example, it can capture a glucuronic acid–like ring.
    uronic_acid_smarts = "C1OC(C(=O)O)C(O)C(O)C1O"
    uronic_query = Chem.MolFromSmarts(uronic_acid_smarts)
    if uronic_query is None:
        return None, None  # Cannot build the query

    # Define heuristic SMARTS for a glycosamine unit.
    # This pattern is for a pyranose ring bearing an amino group (as in glucosamine).
    glycosamine_smarts = "C1OC(C(N)O)C(O)C(O)C1O"
    glyco_query = Chem.MolFromSmarts(glycosamine_smarts)
    if glyco_query is None:
        return None, None

    # Define a SMARTS for a sulfate ester group.
    sulfate_smarts = "OS(=O)(=O)O"  # matches –OS(=O)(=O)O pattern
    sulfate_query = Chem.MolFromSmarts(sulfate_smarts)

    # Find substructure matches.
    uronic_matches = mol.GetSubstructMatches(uronic_query)
    glyco_matches = mol.GetSubstructMatches(glyco_query)
    sulfate_matches = mol.GetSubstructMatches(sulfate_query) if sulfate_query is not None else []

    n_uronic = len(uronic_matches)
    n_glyco = len(glyco_matches)
    n_sulfate = len(sulfate_matches)
    
    # We require at least 2 units of each as the definition involves a polymer.
    if n_uronic < 2 or n_glyco < 2:
        reason = (f"Not enough sugar-like units: found {n_uronic} uronic acid-like "
                  f"and {n_glyco} glycosamine-like rings (need at least 2 of each)")
        return False, reason

    # Check that the counts are roughly equal (difference no more than 1 unit).
    if abs(n_uronic - n_glyco) > 1:
        reason = (f"Counts of repeating units do not alternate properly: "
                  f"{n_uronic} uronic acid-like vs {n_glyco} glycosamine-like rings")
        return False, reason

    # Build classification reason string.
    reason = (f"Found {n_uronic} uronic acid-like and {n_glyco} glycosamine-like rings "
              f"(roughly alternating).")
    if n_sulfate > 0:
        reason += f" Additionally, {n_sulfate} sulfate ester group(s) were detected."

    # Heuristic: if we have the alternating units and optionally sulfate groups, classify as mucopolysaccharide.
    return True, reason