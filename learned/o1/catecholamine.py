"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is defined as 4-(2-aminoethyl)pyrocatechol and derivatives formed by substitution.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for catechol moiety (benzene ring with adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts('c1cc(O)cc(O)c1')
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety found (benzene ring with adjacent hydroxyl groups)"

    # Define SMARTS pattern for aminoethyl side chain attached to benzene ring
    aminoethyl_pattern = Chem.MolFromSmarts('[cH]-[CH2]-[CH2]-[NH2]')
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain attached to ring found"

    # Check if aminoethyl side chain is connected to the catechol ring
    # Find all matches of catechol pattern
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    # Find all matches of aminoethyl pattern
    aminoethyl_matches = mol.GetSubstructMatches(aminoethyl_pattern)
    # Check for connection between catechol and aminoethyl groups
    connected = False
    for cat_match in catechol_matches:
        for amine_match in aminoethyl_matches:
            # The first atom in aminoethyl pattern should be bonded to one of the carbons in catechol ring
            if amine_match[0] in cat_match:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Aminoethyl side chain is not connected to the catechol ring"

    return True, "Contains catechol moiety with aminoethyl side chain attached to ring"