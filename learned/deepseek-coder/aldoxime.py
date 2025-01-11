"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:36586 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is characterized by the presence of the RCH=NOH group, where R is an alkyl or aryl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldoxime pattern: R-CH=NOH, where R is an alkyl or aryl group
    aldoxime_pattern = Chem.MolFromSmarts("[CX3H1](=N[OH])")
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aldoxime group (R-CH=NOH) found"

    # Ensure that the carbon in the R-CH=NOH group is connected to an alkyl or aryl group
    # This helps to avoid false positives where the R group is not appropriate
    for match in mol.GetSubstructMatches(aldoxime_pattern):
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        # Check if the carbon is connected to at least one alkyl or aryl group
        has_alkyl_or_aryl = False
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                if neighbor.GetDegree() >= 1:  # Alkyl or aryl group
                    has_alkyl_or_aryl = True
                    break
        if not has_alkyl_or_aryl:
            return False, "R group in R-CH=NOH is not an alkyl or aryl group"

    return True, "Contains the aldoxime group (R-CH=NOH) with an appropriate R group"