"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is a cholesterol backbone with a fatty acid attached via an ester bond at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more precise SMARTS pattern for the cholesterol core
    # This pattern captures the tetracyclic ring system and the key carbon atoms
    # It's crucial to get the ring fusions right.
    cholesterol_core = Chem.MolFromSmarts("[C]1[C]2[C]([C]([C]3[C]([C]2CC[C]4[C]3[C]CC([C]5[C]4[C])[C]5[C])C)([C])C)[O]")
    
    if cholesterol_core is None:
          return False, "Invalid cholesterol SMARTS"
    # Find matches for the core structure
    matches = mol.GetSubstructMatches(cholesterol_core)

    if not matches:
        return False, "No cholesterol core structure found."
    
    # Check for the ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester group found"

    # Verify that the ester group is attached to the position 3 oxygen of the cholesterol.
    # Find the oxygen from the core and check if its connected to the ester.
    ester_linked_correctly = False
    for core_match in matches:
        oxygen_index_on_core = core_match[-1]
        for ester_match in ester_matches:
            ester_oxygen_index = ester_match[0]
            if mol.GetBondBetweenAtoms(oxygen_index_on_core, ester_oxygen_index) is not None:
                 ester_linked_correctly = True
                 break
        if ester_linked_correctly:
          break
    
    if not ester_linked_correctly:
      return False, "Ester group not attached to position 3 of cholesterol."
    
    # Check for the fatty acid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "No long fatty acid chain attached"
    
    return True, "Cholesterol core with ester group at position 3 and fatty acid chain."