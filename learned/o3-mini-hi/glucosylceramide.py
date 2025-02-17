"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide (cerebrosides in which the monosaccharide head group is glucose)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Glucosylceramides are cerebrosides where the sugar head group is glucose, and the lipid
    portion is bound via an amide bond. Typically the molecule is of high molecular weight 
    (>500 Da).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of the glucose head-group.
    # We use a SMARTS pattern tailored to match a beta-D-glucosyl moiety.
    # Note: Many glucosylceramide structures have a fragment like
    # "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    glucose_smarts = "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    glucose_pattern = Chem.MolFromSmarts(glucose_smarts)
    if glucose_pattern is None:
        return False, "Error in defining glucose SMARTS pattern"
    if not mol.HasSubstructMatch(glucose_pattern):
        # If the exact pattern is not found, try a looser pattern for a glucose ring.
        # This pattern ignores stereochemistry and may pick up a general pyranose ring.
        loose_glucose_smarts = "OC1OC(O)C(O)C(O)C1O"
        loose_glucose_pattern = Chem.MolFromSmarts(loose_glucose_smarts)
        if not mol.HasSubstructMatch(loose_glucose_pattern):
            return False, "Glucose head group not found"

    # Check for the presence of an amide bond,
    # as most ceramide backbones have an amide: C(=O)N.
    amide_smarts = "C(=O)N"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if amide_pattern is None:
        return False, "Error in defining amide SMARTS pattern"
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond detected (expected in the ceramide backbone)"

    # Verify that the molecule is large enough to be a glucosylceramide.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); not typical of a glucosylceramide"

    # Optionally, one might add extra checks such as ensuring the lipid portion has a long carbon chain.
    # For example, count the number of carbon atoms in the non-sugar part of the molecule, or count 
    # rotatable bonds. For our purposes, the sugar and amide checks plus a MW threshold are used.
    
    return True, "Structure has a glucose head group, an amide bond characteristic of a ceramide, and sufficient molecular weight."

# Example test (uncomment below lines to run a test):
# example_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCCC"
# result, reason = is_glucosylceramide(example_smiles)
# print(result, reason)