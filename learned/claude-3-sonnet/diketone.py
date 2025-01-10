"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: diketone
A compound containing exactly two ketone functionalities
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone groups (C=O where C is bonded to two other carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to ensure proper aromatic bond representation
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # SMARTS pattern for a ketone group: carbon double bonded to oxygen,
    # where the carbon is connected to exactly two other carbons
    # [#6] = any carbon
    # X4 = exactly 4 total bonds
    # [CX3](=[OX1])[#6] = carbonyl carbon connected to another carbon
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[#6]")
    
    # Find all ketone matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if len(ketone_matches) != 2:
        return False, f"Found {len(ketone_matches)} ketone groups, need exactly 2"

    # Additional check to exclude cases where we might have matched other carbonyl groups
    # Look for competing functional groups that might have been matched
    
    # Carboxylic acid pattern
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    # Ester pattern
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Amide pattern
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check if any of our ketone matches overlap with these other functional groups
    ketone_carbons = set(match[1] for match in ketone_matches)  # Get the carbonyl carbons
    
    other_carbonyls = set()
    for matches in [acid_matches, ester_matches, amide_matches]:
        for match in matches:
            other_carbonyls.add(match[0])
    
    # If there's overlap between ketone carbons and other carbonyl groups,
    # we might have misidentified some groups
    if ketone_carbons.intersection(other_carbonyls):
        return False, "Contains carbonyl groups that are not ketones"

    # Additional validation - check that the carbons attached to the ketone
    # are not part of other functional groups
    for match in ketone_matches:
        carbonyl_carbon = match[1]
        atom = mol.GetAtomWithIdx(carbonyl_carbon)
        
        # Check neighbors of the carbonyl carbon
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:  # if not carbon
                if neighbor.GetAtomicNum() != 8:  # ignore the oxygen of the ketone
                    return False, "Ketone group is part of another functional group"

    return True, "Contains exactly two ketone groups"