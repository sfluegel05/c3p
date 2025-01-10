"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule belongs to the acrovestone class based on its SMILES string.
    Acrovestone compounds are polyphenols with an isoflavone core structure and various substituents.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule belongs to acrovestone class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more accurate isoflavone core SMARTS pattern
    isoflavone_pattern = Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)c3ccccc23')  # Isoflavone core
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"
    
    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')  # Hydroxyl group pattern
    num_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyl == 0:
        return False, "No hydroxyl groups found"
    
    # Count methoxy groups (-OCH3)
    methoxy_pattern = Chem.MolFromSmarts('CO')  # Methoxy group pattern
    num_methoxy = len(mol.GetSubstructMatches(methoxy_pattern))
    
    # Check for glycosylation (sugar moieties)
    # Using a generic pattern for glycosidic linkage
    glycoside_pattern = Chem.MolFromSmarts('O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O')  # Simple sugar ring
    glycosylated = mol.HasSubstructMatch(glycoside_pattern)
    
    # Additional check for polyphenol character (multiple phenolic hydroxyl groups)
    phenol_pattern = Chem.MolFromSmarts('c[OH]')
    num_phenol = len(mol.GetSubstructMatches(phenol_pattern))
    if num_phenol < 1:
        return False, "No phenolic hydroxyl groups found"
    
    # Acrovestone compounds typically have:
    # - Isoflavone core
    # - Hydroxyl and/or methoxy substituents
    # - Possible glycosylation
    
    reason = "Molecule matches acrovestone class (isoflavone core with appropriate substituents)"
    details = []
    
    details.append(f"Isoflavone core present")
    details.append(f"Number of hydroxyl groups: {num_hydroxyl}")
    if num_methoxy > 0:
        details.append(f"Number of methoxy groups: {num_methoxy}")
    if glycosylated:
        details.append(f"Glycosylation detected")
    else:
        details.append(f"No glycosylation detected")
    details.append(f"Number of phenolic hydroxyl groups: {num_phenol}")
    
    return True, reason + "; " + "; ".join(details)

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'acrovestone',
        'definition': 'A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata can be added as required
}