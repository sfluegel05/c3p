"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are complex polyphenolic compounds with multiple phenolic hydroxyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight (tannins are typically >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Using 400 as lower threshold to catch some smaller tannins
        return False, f"Molecular weight ({mol_wt:.1f}) too low for tannin"

    # Count aromatic rings
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, f"Too few aromatic rings ({aromatic_rings}) for tannin"

    # Look for phenol groups
    phenol_pattern = Chem.MolFromSmarts("[OH]c1ccccc1")
    phenol_matches = len(mol.GetSubstructMatches(phenol_pattern))
    if phenol_matches < 2:
        return False, f"Too few phenol groups ({phenol_matches}) for tannin"

    # Look for galloyl groups (3,4,5-trihydroxybenzoyl)
    galloyl_pattern = Chem.MolFromSmarts("O=C(O)c1c(O)c(O)c(O)cc1")
    catechol_pattern = Chem.MolFromSmarts("Oc1ccc(O)cc1")
    pyrogallol_pattern = Chem.MolFromSmarts("Oc1c(O)c(O)ccc1")
    
    galloyl_matches = len(mol.GetSubstructMatches(galloyl_pattern))
    catechol_matches = len(mol.GetSubstructMatches(catechol_pattern))
    pyrogallol_matches = len(mol.GetSubstructMatches(pyrogallol_pattern))
    
    # Count total hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_count < 3:
        return False, f"Too few hydroxyl groups ({hydroxyl_count}) for tannin"

    # Check for characteristic structural features
    has_galloyl = galloyl_matches > 0
    has_catechol = catechol_matches > 0
    has_pyrogallol = pyrogallol_matches > 0
    
    if not (has_galloyl or has_catechol or has_pyrogallol):
        return False, "Missing characteristic tannin structural features (galloyl/catechol/pyrogallol groups)"

    # Look for ester bonds (common in hydrolyzable tannins)
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)O-[#6]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    
    # Calculate ring complexity
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Classify based on combined features
    if (hydroxyl_count >= 3 and 
        (has_galloyl or (has_catechol and has_pyrogallol)) and
        ring_count >= 2 and
        mol_wt >= 400):
        
        features = []
        if has_galloyl:
            features.append(f"{galloyl_matches} galloyl groups")
        if has_catechol:
            features.append(f"{catechol_matches} catechol groups")
        if has_pyrogallol:
            features.append(f"{pyrogallol_matches} pyrogallol groups")
        if ester_matches > 0:
            features.append(f"{ester_matches} ester bonds")
            
        return True, f"Contains {hydroxyl_count} hydroxyl groups, {ring_count} rings, and {', '.join(features)}"

    return False, "Does not meet minimum structural requirements for tannin classification"