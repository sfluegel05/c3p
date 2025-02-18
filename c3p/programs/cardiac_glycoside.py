"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: Cardiac glycoside 
Definition: Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.
This program attempts to determine if the molecule (given by SMILES) has:
  1) A steroid nucleus (fused cyclohexane/cyclopentane system)
  2) A lactone (butenolide) ring system
  3) At least one sugar residue (e.g. a pyranose ring pattern)
Note: SMARTS definition is approximate, so false negatives/positives are possible.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Cardiac glycosides are defined as steroid lactones containing sugar residues.
    This function checks for:
      - A steroid nucleus (approximate fused tetracyclic system)
      - A lactone ring (furan-2-one like substructure)
      - At least one sugar (cyclic acetal with multiple hydroxyls) attached
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a cardiac glycoside, otherwise False.
        str: Reason for classification outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a steroid nucleus.
    # This is a rough approximation for the cyclopentanoperhydrophenanthrene system.
    steroid_smarts = "[$([C@@H]1CC[C@H]2[C@@H]3CC[C@]4(C)[C@H](CC[C@]4(C)C3)CC[C@]12C)]"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus (cyclopentanoperhydrophenanthrene system) found"
    
    # Define a SMARTS pattern for a butenolide lactone (furan-2-one type)
    # This pattern looks for a 5-membered ring containing a carbonyl and an oxygen.
    lactone_smarts = "O=C1C=CCO1"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group (butenolide/furan-2-one) found"

    # Define a SMARTS pattern for a sugar residue. We use an approximate hexopyranose pattern.
    sugar_smarts = "[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar residue (pyranose-like ring) found"

    # Optionally, one can check for the molecular size expected for cardiac glycosides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low to be a typical cardiac glycoside"

    return True, "Contains a steroid nucleus, a lactone (butenolide) ring and at least one sugar residue"