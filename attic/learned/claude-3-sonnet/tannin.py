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
        # Core phenolic patterns (more general)
        'phenol': ['[OH]c1ccccc1', '[OH]c1cccc([OH])c1', '[OH]c1cc([OH])ccc1'],
        'catechol': ['[OH]c1c([OH])cccc1', '[OH]c1c([OH])cc([OH])cc1'],
        'pyrogallol': ['[OH]c1c([OH])c([OH])ccc1', '[OH]c1c([OH])c([OH])cc([OH])c1'],
        
        # Ester and ether linkages (more flexible)
        'ester': ['C(=O)O[CH0,CH1,CH2]', 'O=C([O,N])c'],
        'ether': ['cOc', 'COc', 'COC'],
        
        # Sugar-related patterns (more general)
        'sugar_core': ['[OH]C1O[CH]([CH]([OH])[CH]([OH])[CH]1[OH])',
                      'OC1OC([CH]([OH])[CH]([OH]))CC1',
                      'C1OC([CH]O)C([OH])C([OH])C1'],
        
        # General features
        'hydroxyl': ['[OH]'],
        'aromatic': ['c1ccccc1'],
        'carbonyl': ['C(=O)'],
        
        # Complex linkage patterns
        'flavonoid_core': ['O1c2c(cc([OH])cc2)CC(c2ccc([OH])cc2)C1',
                          'Oc1cc(O)c2C(=O)CC(c3ccc(O)cc3)Oc2c1'],
        
        # Oligomeric patterns
        'biaryl': ['c1ccc(-c2ccccc2)cc1',
                   'c1cc(-c2cc([OH])c([OH])cc2)ccc1']
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
    if counts['hydroxyl'] < 3:
        return False, f"Too few hydroxyl groups ({counts['hydroxyl']}) for tannin"

    ring_count = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    
    if ring_count < 2:
        return False, f"Too few rings ({ring_count}) for tannin"
    if aromatic_rings < 1:
        return False, "No aromatic rings found"

    # Calculate phenolic character
    phenolic_groups = counts['phenol'] + counts['catechol'] + counts['pyrogallol']
    if phenolic_groups < 2:
        return False, f"Insufficient phenolic groups ({phenolic_groups})"

    # Verify polyphenolic nature
    num_carbons = len([a for a in mol.GetAtoms() if a.GetAtomicNum() == 6])
    phenolic_density = phenolic_groups / num_carbons
    if phenolic_density < 0.05:  # Lowered threshold to 5%
        return False, "Insufficient phenolic group density"

    # Check for hydrolyzable tannin characteristics
    is_hydrolyzable = (counts['ester'] >= 1 and 
                      (counts['sugar_core'] >= 1 or 
                       (counts['carbonyl'] >= 2 and counts['ether'] >= 1)))

    # Check for condensed tannin characteristics
    is_condensed = ((counts['flavonoid_core'] >= 1) or
                   (counts['biaryl'] >= 1 and phenolic_groups >= 3) or
                   (counts['catechol'] + counts['pyrogallol'] >= 3 and counts['ether'] >= 1))

    # Build feature description
    features = []
    if counts['phenol'] > 0:
        features.append(f"{counts['phenol']} phenolic groups")
    if counts['catechol'] > 0:
        features.append(f"{counts['catechol']} catechol groups")
    if counts['pyrogallol'] > 0:
        features.append(f"{counts['pyrogallol']} pyrogallol groups")
    if counts['ester'] > 0:
        features.append(f"{counts['ester']} ester bonds")
    if counts['sugar_core'] > 0:
        features.append("sugar core")
    if counts['flavonoid_core'] > 0:
        features.append("flavonoid core")

    # Final classification
    if is_hydrolyzable or is_condensed:
        tannin_type = []
        if is_hydrolyzable:
            tannin_type.append("hydrolyzable")
        if is_condensed:
            tannin_type.append("condensed")
        
        return True, f"Contains {counts['hydroxyl']} hydroxyl groups, {ring_count} rings, and {', '.join(features)}. Classified as {' and '.join(tannin_type)} tannin"

    # Additional check for complex polyphenols that might not fit standard patterns
    if (phenolic_groups >= 4 and 
        ring_count >= 3 and 
        (counts['ether'] >= 2 or counts['ester'] >= 2)):
        return True, f"Complex polyphenolic structure with {counts['hydroxyl']} hydroxyl groups and {ring_count} rings"

    return False, "Does not meet structural requirements for tannin classification"