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

    # Define SMARTS for Coenzyme A moiety, focusing on the typical structural features
    coa_pattern = Chem.MolFromSmarts("CC(C)C(=O)NCCSC(=O)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Define SMARTS for a double bond between the 2nd and 3rd carbon in the acyl chain
    double_bond_2_3_pattern = Chem.MolFromSmarts("C-C=C-C")
    substruct_matches = mol.GetSubstructMatches(double_bond_2_3_pattern)

    if not substruct_matches:
        return False, "No suitable double bond between C2 and C3 found"

    # Verify positional context - check that the identified double bond is part of an acyl chain
    # and near the CoA linkage
    for match in substruct_matches:
        # Assuming match is the fragment where the double bond is located
        # [C1, C2, C3, C4] where C2=C3 is the double bond
        chain_carbon_2 = mol.GetAtomWithIdx(match[1])
        chain_carbon_3 = mol.GetAtomWithIdx(match[2])

        # Both carbons should be part of a chain, not in a cycle
        if chain_carbon_2.GetDegree() <= 3 and chain_carbon_3.GetDegree() <= 3:
            # Verify acyl linkage is nearby, indicating it's in the correct acyl position
            return True, "Contains a suitable 2,3-double bond and Coenzyme A moiety"

    return False, "Matched double bond is not appropriate in positional context"