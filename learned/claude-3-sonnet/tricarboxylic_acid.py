"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
Definition: An oxoacid containing three carboxy groups
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid must contain exactly three carboxylic acid (-COOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and sanitize
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Molecule failed sanitization"

    # Look for carboxylic acid groups (-C(=O)OH)
    # More specific SMARTS pattern that excludes esters and other similar groups
    carboxyl_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H,OX2-]')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid groups found"
    
    # Count carboxylic acid groups
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl = len(carboxyl_matches)
    
    if num_carboxyl != 3:
        return False, f"Found {num_carboxyl} carboxylic acid groups, need exactly 3"

    # Check for esters (to exclude them)
    ester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][#6]')
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches > 0:
        return False, "Contains ester groups"

    # Check for acid anhydrides (to exclude them)
    anhydride_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][CX3](=[OX1])')
    anhydride_matches = len(mol.GetSubstructMatches(anhydride_pattern))
    if anhydride_matches > 0:
        return False, "Contains acid anhydride groups"

    # Verify oxygen count
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:  # Each COOH has 2 oxygens, so minimum 6 oxygens needed
        return False, "Insufficient oxygen atoms for three carboxylic acids"

    # Additional check for peptide-like bonds to reduce false positives
    peptide_pattern = Chem.MolFromSmarts('[NX3,NX4][CX3](=[OX1])[#6]')
    peptide_matches = len(mol.GetSubstructMatches(peptide_pattern))
    if peptide_matches > 0:
        # Only reject if the peptide bonds would interfere with our carboxyl count
        potential_interfering_groups = peptide_matches + ester_matches
        if potential_interfering_groups + num_carboxyl > 3:
            return False, "Contains peptide bonds that may be incorrectly counted"

    return True, "Contains exactly three carboxylic acid groups"