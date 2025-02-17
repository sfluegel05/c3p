"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide 
Definition: any cerebroside where the monosaccharide head group is glucose.
We require that the molecule contains a single glucose ring (attached via an oxygen),
an amide bond that connects the fatty acid to the sphingoid base, and a molecular 
weight typical of a ceramide (>500 Da).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    
    The classifier checks for:
      1) A glucose head-group (using a SMARTS pattern that is both specific and loose).
         We also require that exactly one such sugar ring is found.
      2) An amide bond (using a somewhat general pattern "[CX3](=O)[NX3]").
      3) A molecular weight >500 Da.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is identified as a glucosylceramide, False otherwise.
      str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Check for a single glucose head-group -----
    # A strict pattern that often appears in beta-D-glucosyl groups attached via an oxygen.
    strict_glucose_smarts = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O"
    strict_glucose_pattern = Chem.MolFromSmarts(strict_glucose_smarts)
    if strict_glucose_pattern is None:
        return False, "Error in defining strict glucose SMARTS"
    glucose_matches = mol.GetSubstructMatches(strict_glucose_pattern)
    
    # If not found using the strict pattern, try a looser one (ignoring stereochemistry).
    if not glucose_matches:
        loose_glucose_smarts = "OC1OC(O)C(O)C(O)C1O"
        loose_glucose_pattern = Chem.MolFromSmarts(loose_glucose_smarts)
        if loose_glucose_pattern is None:
            return False, "Error in defining loose glucose SMARTS"
        glucose_matches = mol.GetSubstructMatches(loose_glucose_pattern)
    
    num_glucose = len(glucose_matches)
    if num_glucose == 0:
        return False, "Glucose head group not found"
    if num_glucose > 1:
        return False, f"Found {num_glucose} sugar rings; expected exactly one glucose head group"
    
    # ----- Check for the presence of an amide bond -----
    # Use a more general pattern that catches most amide bonds.
    amide_smarts = "[CX3](=O)[NX3]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if amide_pattern is None:
        return False, "Error in defining amide SMARTS pattern"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 1:
        return False, "No amide bond detected (expected in the ceramide backbone)"
    
    # ----- Check molecular weight (typical glucosylceramides are large) -----
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); not typical of a glucosylceramide"
    
    # If all criteria are met, return True.
    return True, "Structure has a single glucose head group, an amide bond characteristic of a ceramide, and sufficient molecular weight."

# Example test (uncomment to run):
# example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCCC"
# result, reason = is_glucosylceramide(example_smiles)
# print(result, reason)