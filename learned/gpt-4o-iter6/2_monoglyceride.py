"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is defined by an ester linkage at the second carbon of a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-monoglyceride pattern:
    # The substructure is described as:
    # C1 is the first carbon with a hydroxyl group, C2 is the second carbon with 
    # the ester linkage, and C3 is the third carbon with another hydroxyl group
    glycerol_2_mono_pattern = Chem.MolFromSmarts("C(CO)C(OC(=O))CO")
    
    # Validate that the ester is specifically linked to the second carbon
    matches = mol.GetSubstructMatches(glycerol_2_mono_pattern)
    if not matches:
        return False, "No specific 2-monoglyceride pattern found with esterification at C2"

    for match in matches:
        c2_index = match[1]  # Index of C2 in the matches
        c1, c2, c3 = mol.GetAtomWithIdx(c2_index).GetNeighbors()
        # Confirm C1 and C3 are the carbon atoms with hydroxyl groups, not oxygens
        if (c1.GetSymbol() == 'C' and c3.GetSymbol() == 'C' and 
            any(n.GetSymbol() == 'O' for n in c1.GetNeighbors()) and 
            any(n.GetSymbol() == 'O' for n in c3.GetNeighbors())):
            return True, "Contains 2-monoglyceride structure with acyl group esterified at the second position"
    
    return False, "The structure does not match the specific pattern of a 2-monoglyceride"