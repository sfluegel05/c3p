"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a glucuronic acid moiety with a glycosidic bond at the C1 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glucuronic acid substructure using SMARTS (less strict stereochemistry)
    # This pattern checks for a six-membered ring, with an oxygen in the ring,
    # a carboxyl group attached to one carbon and two hydroxyl groups
    # It does not check for stereochemistry explicitly
    glucuronic_acid_smarts = "C1OC(C(O)C(O)C1C(=O)O)"
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)


    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
         return False, "No glucuronic acid substructure found"

    # Get the matches for the substructure
    matches = mol.GetSubstructMatches(glucuronic_acid_pattern)

    #Check that one of the C1 is attached to an oxygen
    glycosidic_bond_found = False
    for match in matches:
       c1_index = match[0] # the index of C1 in the molecule according to the SMARTS
       c1_atom = mol.GetAtomWithIdx(c1_index)
       for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in [a.GetIdx() for a in c1_atom.GetNeighbors() if a.GetAtomicNum() == 6]: # check for oxygen NOT part of the ring
               glycosidic_bond_found = True
               break
       if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
       return False, "No glycosidic bond at C1 of glucuronic acid."

    return True, "Contains beta-D-glucuronic acid with a glycosidic bond at the C1 position."