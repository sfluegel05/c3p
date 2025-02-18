"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo fatty acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group 
            in addition to the carboxylic acid group.
Additional criteria:
  - The molecule must have exactly one carboxylic acid group.
  - The molecule must be acyclic (no rings) and have a sufficient number of carbons.
  - One extra oxo (ketone or aldehyde) group must be present outside that acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    The molecule must contain a carboxylic acid group (as present in fatty acids)
    and an additional ketone or aldehyde function elsewhere.
    Furthermore, the molecule is expected to be fatty-acid like,
    meaning it should be acyclic, contain only one acid group, and have a minimum number of carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oxo fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion 1: Check for carboxylic acid group.
    # SMARTS for a carboxylic acid: a carbonyl carbon (sp2) double-bonded to O and single-bonded to an -OH.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) == 0:
        return False, "Molecule does not contain a carboxylic acid group."
    elif len(acid_matches) > 1:
        # Many fatty acids are monocarboxylic. More than one acid group suggests it is not a typical fatty acid.
        return False, "Molecule contains more than one carboxylic acid group, not a typical fatty acid."
    
    # Save the carbonyl carbon index of the acid group.
    acid_carbonyl_indices = {match[0] for match in acid_matches}
    
    # Criterion 2: Check that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and does not appear to be a linear fatty acid."

    # Criterion 3: Check for a minimum carbon count (at least 5 carbons).
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, f"Molecule contains too few carbons ({carbon_count}) to be a fatty acid."

    # Criterion 4: Look for an extra oxo group in the form of a ketone or aldehyde (outside the acid group).
    extra_oxo_found = False

    # Check for ketone group.
    # SMARTS for ketone: carbonyl carbon bonded to two carbon atoms.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    for match in mol.GetSubstructMatches(ketone_pattern):
        # match is a tuple (R, C, R') where C is the carbonyl carbon.
        carbonyl_index = match[1]
        if carbonyl_index not in acid_carbonyl_indices:
            extra_oxo_found = True
            break

    # If not found as a ketone, try finding an aldehyde group.
    if not extra_oxo_found:
        # SMARTS for an aldehyde: the carbonyl carbon is bonded to one carbon and one hydrogen.
        aldehyde_smarts = "[#6][CX3H1](=O)"
        aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            carbonyl_index = match[1]
            if carbonyl_index not in acid_carbonyl_indices:
                extra_oxo_found = True
                break

    if not extra_oxo_found:
        return False, "Molecule does not contain an additional aldehydic or ketonic group outside the carboxylic acid."

    return True, "Molecule is an oxo fatty acid: it is acyclic, has a single carboxylic acid group and an extra oxo function."