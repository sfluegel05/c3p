"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo fatty acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group 
            in addition to the carboxylic acid group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    The molecule must contain a carboxylic acid group (as found in fatty acids)
    and an additional carbonyl moiety in the form of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 1. Check for the presence of a carboxylic acid group.
    # The SMARTS below looks for a carbon atom (C, sp2 hybridized) double-bonded to oxygen
    # and singly bonded to an -OH group. More robust patterns could be used but this is sufficient.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Molecule does not contain a carboxylic acid group."

    # Retrieve indices of the carboxylic acid carbonyl atom(s)
    # (we assume the first atom in the match is the carbonyl carbon)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    acid_carbonyl_indices = {match[0] for match in acid_matches}
    
    # 2. Check for an additional oxo group in the form of a ketone or aldehyde.
    # Define SMARTS for ketone: [#6][CX3](=O)[#6] matches a carbonyl group flanked by carbons.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    extra_oxo_found = False
    for match in ketone_matches:
        # match is a tuple of three atoms: (R-carbon, carbonyl carbon, R'-carbon)
        carbonyl_index = match[1]
        # Ensure that this carbonyl is not the one in the carboxylic acid group.
        if carbonyl_index not in acid_carbonyl_indices:
            extra_oxo_found = True
            break
    
    # If not found as a ketone, check if an aldehyde group exists.
    if not extra_oxo_found:
        # Define SMARTS for aldehyde: [#6][CX3H1](=O) - one alkyl group attached 
        # and one hydrogen present on the carbonyl center.
        aldehyde_smarts = "[#6][CX3H1](=O)"
        aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
        for match in aldehyde_matches:
            carbonyl_index = match[1]
            if carbonyl_index not in acid_carbonyl_indices:
                extra_oxo_found = True
                break

    # If neither an extra ketone nor an extra aldehyde group was found, it is not an oxo fatty acid.
    if not extra_oxo_found:
        return False, "Molecule does not contain an additional aldehydic or ketonic group outside the carboxylic acid."
    
    # Passed both criteria: carboxylic acid and additional oxo function.
    return True, "Molecule contains a carboxylic acid and an extra oxo (ketone or aldehyde) group, classifying it as an oxo fatty acid."