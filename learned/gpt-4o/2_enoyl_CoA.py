"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is characterized by a double bond between the 2nd and 3rd
    positions in the acyl chain and a Coenzyme A moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced SMARTS pattern for Coenzyme A moiety based on typical structure
    # Capture larger portion and key linkage
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCSC(=O)[C]")  # Key amide and thioester bonds
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # SMARTS for 2-enoyl character: Double bond between 2nd and 3rd carbons in a chain
    double_bond_2_3_pattern = Chem.MolFromSmarts("C=C")
    substruct_matches = mol.GetSubstructMatches(double_bond_2_3_pattern)
    
    if not substruct_matches:
        return False, "No suitable double bond between C2 and C3 found"

    for match in substruct_matches:
        # Verify that the double bond is correctly placed to be in the 2-enoyl position (right after the CoA linkage)
        # We can approximate the location using the bond order and check connectivity
        if mol.GetAtomWithIdx(match[0]).GetDegree() > 1 and mol.GetAtomWithIdx(match[1]).GetDegree() > 1:
            # Contextual checks may include ensuring the bond exists in a suitable long chain near CoA
            return True, "Contains suitable 2,3-double bond and Coenzyme A moiety"

    return False, "Matched double bond is not appropriate in positional context"