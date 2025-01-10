"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is defined by an ester linkage at the second carbon of a glycerol backbone,
    leaving the first and third carbons with hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a 2-monoglyceride
    # The pattern looks for the glycerol core with an ester linkage at the second carbon:
    # [OH]-C([OH])COC(=O)
    # This pattern checks for the ester at C2 specifically with C1 and C3 having hydroxyls
    glycerol_2_mono_pattern = Chem.MolFromSmarts("[OH]C(COC(=O))CO")

    # Search for matches
    matches = mol.GetSubstructMatches(glycerol_2_mono_pattern)
    
    if not matches:
        return False, "No 2-monoglyceride pattern found with esterification at C2"

    # Confirm the ester linkage is correctly situated
    for match in matches:
        c1_index, c2_index, ester_oxygen_index = match[:3]

        # Check bonding for C1 (should have hydroxyl group)
        c1 = mol.GetAtomWithIdx(c1_index)
        if not any(neigh.GetAtomicNum() == 8 for neigh in c1.GetNeighbors()):
            return False, "C1 does not have a hydroxyl group"

        # Check bonding for C2 (should be linked to the ester oxygen)
        ester_o = mol.GetAtomWithIdx(ester_oxygen_index)
        if not any(neigh.GetIdx() == c2_index for neigh in ester_o.GetNeighbors()):
            return False, "Ester linkage not at C2"

    return True, "Contains the structure of 2-monoglyceride with ester at C2"