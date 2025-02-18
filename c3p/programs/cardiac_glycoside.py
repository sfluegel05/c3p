"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside 
Definition: Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.
This program checks for the presence of:
  1) A steroid nucleus (approximated as a fused tetracyclic skeleton)
  2) A lactone ring (approximated as a butenolide, furan-2-one substructure)
  3) At least one sugar residue (detected via a pyranose-like ring; if a strict pattern fails, a looser one is attempted)
Note: The SMARTS definitions are approximate and may not capture every nuance of the structural class.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.

    Cardiac glycosides are steroid lactones with sugar moieties.
    This function checks for:
      - A steroid nucleus (fused tetracyclic system, roughly representing the cyclopentanoperhydrophenanthrene system)
      - A lactone (butenolide/furan-2-one) ring system
      - At least one sugar residue (pyranose-like ring)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a cardiac glycoside, False otherwise.
        str: Reason for the classification outcome.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Check for the steroid nucleus using a simplified SMARTS pattern that represents a fused tetracyclic system.
    # This pattern is a rough approximation: a set of three six-membered rings and one five-membered ring fused together.
    steroid_smarts = "C1CC2CC3CC(C1)C2CC3"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus (fused tetracyclic system) found"

    # (2) Look for a lactone group.
    # The SMARTS here targets a butenolide (furan-2-one) ring: a 5-membered ring containing a carbonyl and a C=C unit.
    lactone_smarts = "O=C1C=CCO1"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group (butenolide/furan-2-one) found"

    # (3) Look for a sugar residue.
    # First try a stereochemically defined pyranose-like ring.
    sugar_smarts = "[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if not mol.HasSubstructMatch(sugar_pattern):
        # Fallback: use a looser pattern for a six-membered ring containing two oxygens.
        sugar_smarts_loose = "O1CC(O)CC(O)C1"
        sugar_pattern_loose = Chem.MolFromSmarts(sugar_smarts_loose)
        if not mol.HasSubstructMatch(sugar_pattern_loose):
            return False, "No sugar residue (pyranose-like ring) found"

    # Optionally check that the molecule is not too small.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for typical cardiac glycosides"

    return True, "Contains a steroid nucleus, a lactone ring, and at least one sugar residue"