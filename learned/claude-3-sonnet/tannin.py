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

    # Check molecular weight (tannins are typically >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for tannin"

    # Count aromatic rings
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, f"Too few aromatic rings ({aromatic_rings}) for tannin"

    # Define structural patterns
    patterns = {
        'galloyl': ['O=C(O)c1c(O)c(O)c(O)cc1', 'O=C([O,N])c1c(O)c(O)c(O)cc1'],
        'catechol': ['Oc1ccc(O)cc1', 'Oc1cccc(O)c1'],
        'pyrogallol': ['Oc1c(O)c(O)ccc1', 'Oc1cc(O)c(O)cc1'],
        'flavonoid': ['O=C1CC(c2ccccc2)Oc2ccccc12'],
        'ester': ['C(=O)O[CH0]', 'C(=O)OC'],
        'hydroxyl': ['[OH]'],
        'linked_flavonoid': ['O1c2c(cc(O)cc2)C(c2ccc(O)cc2)C[C@H]1c1c(O)cc(O)c2c1O[C@H](c1ccc(O)cc1)CC2']
    }
    
    # Count structural features
    counts = {}
    for name, smarts_list in patterns.items():
        counts[name] = 0
        for smarts in smarts_list:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                counts[name] += len(mol.GetSubstructMatches(pattern))
    
    # Basic requirements
    if counts['hydroxyl'] < 5:
        return False, f"Too few hydroxyl groups ({counts['hydroxyl']}) for tannin"

    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3:
        return False, f"Too few rings ({ring_count}) for tannin"

    # Check for characteristic tannin features
    is_hydrolyzable = (counts['galloyl'] > 0 and counts['ester'] > 0)
    is_condensed = (counts['catechol'] + counts['pyrogallol'] >= 2) or counts['linked_flavonoid'] > 0
    has_phenolic_groups = (counts['catechol'] + counts['pyrogallol'] + counts['galloyl'] >= 2)

    if not (is_hydrolyzable or is_condensed) or not has_phenolic_groups:
        return False, "Insufficient tannin structural features"

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
    if counts['linked_flavonoid'] > 0:
        features.append("linked flavonoid units")

    # Final classification
    if ((is_hydrolyzable or is_condensed) and 
        has_phenolic_groups and 
        counts['hydroxyl'] >= 5 and 
        ring_count >= 3):
        
        return True, f"Contains {counts['hydroxyl']} hydroxyl groups, {ring_count} rings, and {', '.join(features)}"

    return False, "Does not meet minimum structural requirements for tannin classification"