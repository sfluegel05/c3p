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

    # Define a more comprehensive SMARTS pattern for Coenzyme A
    # This pattern captures a broader section of the CoA moiety,
    # with focus on the adenosine part attached via a phosphothioester linkage
    coenzyme_a_pattern = Chem.MolFromSmarts("C=O)SCCNC(=O)CC")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "No Coenzyme A moiety found"

    # Define a pattern for a flexible double bond starting early in the acyl chain,
    # but not assuming it's directly after CoA, to handle branched structures
    double_bond_2_3_pattern = Chem.MolFromSmarts("C=CC")
    matches = mol.GetSubstructMatches(double_bond_2_3_pattern)

    if not matches:
        return False, "No double bond between the 2nd and 3rd positions"

    # Validate positional match within the acyl chain near the CoA linkage
    # Ensure the double bond isn't part of a cyclized structure, irrelevant to the required connectivity
    for match in matches:
        # Check the atom indices of the double bond match to ensure it's within the acyl chain directly bonded to CoA
        start, end = match[0], match[1]
        atom_start = mol.GetAtomWithIdx(start)
        atom_end = mol.GetAtomWithIdx(end)

        # Ensure it's a linear acyl chain fragment for the CoA moiety
        # This is a simple heuristic check; more complex molecular graph analysis would be better in a production system.
        if atom_start.GetDegree() <= 2 and atom_end.GetDegree() <= 2:
            return True, "Contains a 2,3-double bond and Coenzyme A moiety"

    return False, "Matched double bond is not appropriate in positional context"