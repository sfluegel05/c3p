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

    # More specific SMARTS patterns for carboxylic acids
    # Includes both protonated and deprotonated forms
    carboxyl_patterns = [
        '[CX3](=[OX1])[OX2H]',  # Standard COOH
        '[CX3](=[OX1])[OX2-]',  # Deprotonated form
        '[CX3](=[OX1])[OX2]([H])',  # Explicit H form
    ]
    
    total_carboxyl_count = 0
    for pattern in carboxyl_patterns:
        pat = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pat)
        total_carboxyl_count += len(matches)
    
    if total_carboxyl_count == 0:
        return False, "No carboxylic acid groups found"
    
    if total_carboxyl_count != 3:
        return False, f"Found {total_carboxyl_count} carboxylic acid groups, need exactly 3"

    # Check for potentially interfering groups
    ester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][CX4]')  # More specific ester pattern
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    
    # Refined anhydride pattern
    anhydride_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][CX3](=[OX1])')
    anhydride_matches = len(mol.GetSubstructMatches(anhydride_pattern))
    
    # Check for acid chlorides
    acid_chloride_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[Cl]')
    acid_chloride_matches = len(mol.GetSubstructMatches(acid_chloride_pattern))
    
    # More specific peptide pattern
    peptide_pattern = Chem.MolFromSmarts('[NX3][CX3](=[OX1])[#6]')
    peptide_matches = len(mol.GetSubstructMatches(peptide_pattern))

    # Calculate total oxygen atoms in carboxylic acid groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Various rejection conditions
    if ester_matches > 0 and ester_matches + total_carboxyl_count > 3:
        return False, "Contains ester groups that interfere with carboxylic acid count"
    
    if anhydride_matches > 0:
        return False, "Contains acid anhydride groups"
        
    if acid_chloride_matches > 0:
        return False, "Contains acid chloride groups"
    
    if o_count < 6:
        return False, "Insufficient oxygen atoms for three carboxylic acids"

    # Only reject peptide-containing structures if they would lead to miscounting
    if peptide_matches > 0:
        # Count carbonyls to ensure we're not double-counting
        carbonyl_pattern = Chem.MolFromSmarts('[CX3]=[OX1]')
        carbonyl_count = len(mol.GetSubstructMatches(carbonyl_pattern))
        if carbonyl_count > total_carboxyl_count + peptide_matches:
            return False, "Contains peptide bonds that interfere with carboxylic acid count"

    return True, "Contains exactly three carboxylic acid groups"