"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are complex polyphenolic compounds including both hydrolyzable and condensed types.
    
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

    # Check molecular weight (tannins typically >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for tannin"

    # Define structural patterns
    patterns = {
        # Core structural patterns
        'galloyl': ['O=C(O)c1c(O)c(O)c(O)cc1', 'O=C([O,N])c1c(O)c(O)c(O)cc1'],
        'catechol': ['[OH]c1[cX3]cc([OH])cc1', '[OH]c1[cX3]cccc1[OH]'],
        'pyrogallol': ['[OH]c1c([OH])c([OH])[cX3]cc1', '[OH]c1[cX3]c([OH])c([OH])cc1'],
        'ellagic_core': ['O=C1Oc2c(O)c(O)cc3c2c2c(c1)c(O)c(O)cc2oc13'],
        'glucose_core': ['[OH]C1OC(CO)[CH]([OH])[CH]([OH])[CH]1[OH]'],
        
        # Linkage patterns
        'ester': ['C(=O)O[CH0]', 'C(=O)OC'],
        'flavonoid_link': ['O1c2c(cc([OH])cc2)C(c2ccc([OH])cc2)CC1'],
        'c4c8_link': ['c1c(c2)c(O)cc(c1)C(O)Cc1c(O)cc(O)c(c1)C2'],
        
        # General features
        'hydroxyl': ['[OH]'],
        'aromatic_oxygen': ['c1ccc(Oc2)cc1', 'c1cc(Oc2)ccc1']
    }
    
    # Count structural features
    counts = {}
    matches = {}
    for name, smarts_list in patterns.items():
        counts[name] = 0
        matches[name] = []
        for smarts in smarts_list:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                these_matches = mol.GetSubstructMatches(pattern)
                counts[name] += len(these_matches)
                matches[name].extend(these_matches)

    # Basic requirements
    if counts['hydroxyl'] < 4:
        return False, f"Too few hydroxyl groups ({counts['hydroxyl']}) for tannin"

    ring_count = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    
    if ring_count < 2:
        return False, f"Too few rings ({ring_count}) for tannin"
    if aromatic_rings < 1:
        return False, "No aromatic rings found"

    # Check for characteristic tannin features
    is_hydrolyzable = False
    is_condensed = False
    
    # Hydrolyzable tannin check
    if counts['galloyl'] > 0 and counts['ester'] > 0:
        if counts['glucose_core'] > 0 or counts['ellagic_core'] > 0:
            is_hydrolyzable = True
            
    # Condensed tannin check
    if (counts['flavonoid_link'] > 0 or counts['c4c8_link'] > 0):
        if (counts['catechol'] + counts['pyrogallol'] >= 2):
            is_condensed = True
    elif counts['catechol'] + counts['pyrogallol'] >= 3:
        # Check if catechol/pyrogallol groups are properly connected
        if counts['aromatic_oxygen'] >= 1:
            is_condensed = True

    # Verify phenolic group density
    num_carbons = len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 6])
    phenolic_density = (counts['catechol'] + counts['pyrogallol'] + counts['galloyl']) / num_carbons
    if phenolic_density < 0.1:  # At least 10% of carbons should be part of phenolic groups
        return False, "Insufficient phenolic group density"

    # Build description of features
    features = []
    if counts['galloyl'] > 0:
        features.append(f"{counts['galloyl']} galloyl groups")
    if counts['catechol'] > 0:
        features.append(f"{counts['catechol']} catechol groups")
    if counts['pyrogallol'] > 0:
        features.append(f"{counts['pyrogallol']} pyrogallol groups")
    if counts['ester'] > 0:
        features.append(f"{counts['ester']} ester bonds")
    if counts['ellagic_core'] > 0:
        features.append("ellagic acid core")
    if counts['glucose_core'] > 0:
        features.append("glucose core")

    # Final classification
    if is_hydrolyzable or is_condensed:
        tannin_type = []
        if is_hydrolyzable:
            tannin_type.append("hydrolyzable")
        if is_condensed:
            tannin_type.append("condensed")
        
        return True, f"Contains {counts['hydroxyl']} hydroxyl groups, {ring_count} rings, and {', '.join(features)}. Classified as {' and '.join(tannin_type)} tannin"

    return False, "Does not meet structural requirements for tannin classification"