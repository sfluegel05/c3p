"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: polymer
A mixture composed of macromolecules of different kinds which may be differentiated by composition, 
length, degree of branching etc.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and check for multiple components
    components = smiles.split('.')
    if len(components) == 1:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
        mols = [mol]
    else:
        mols = [Chem.MolFromSmiles(c) for c in components]
        if any(m is None for m in mols):
            return False, "Invalid SMILES string"

    # Analyze molecular properties of components
    component_props = []
    total_mw = 0
    for mol in mols:
        props = {
            'mw': Descriptors.ExactMolWt(mol),
            'atoms': mol.GetNumAtoms(),
            'rings': rdMolDescriptors.CalcNumRings(mol),
            'rot_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'charge': Chem.GetFormalCharge(mol)
        }
        component_props.append(props)
        total_mw += props['mw']

    # Check for polymer characteristics
    num_components = len(mols)
    
    # Look for common polymer patterns
    polymer_patterns = [
        ('[*]-[*]-[*]-[*]-[*]-[*]', 'linear chain'),
        ('C(=O)O', 'carboxylic acid'),
        ('C(=O)N', 'amide'),
        ('COC', 'ether linkage'),
        ('CC(C)(C)C', 'branched alkyl')
    ]
    
    pattern_counts = {}
    for pattern, name in polymer_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            pattern_counts[name] = sum(len(mol.GetSubstructMatches(patt)) for mol in mols)

    # Classification logic
    if num_components >= 2:
        # Check for salt forms and charged species
        charges = [p['charge'] for p in component_props]
        if any(c != 0 for c in charges):
            large_components = sum(1 for p in component_props if p['mw'] > 300)
            if large_components >= 1:
                return True, "Multi-component ionic polymer system"
        
        # Check for repeating structural motifs
        if any(count >= 3 for count in pattern_counts.values()):
            return True, "Contains repeating polymer subunits"
        
        # Check for size diversity in components
        mw_std = sum((p['mw'] - total_mw/num_components)**2 for p in component_props)**0.5
        if mw_std > 100 and num_components >= 3:
            return True, "Mixture of diverse macromolecular components"

    # Single component analysis
    for mol, props in zip(mols, component_props):
        # Check for very large molecules with repeating units
        if props['mw'] > 800 and props['rot_bonds'] > 15:
            # Look for repeating units
            for patt_name, count in pattern_counts.items():
                if count >= 4:
                    return True, f"Large molecule with repeating {patt_name} units"
        
        # Check for branched structures
        branch_points = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[*]([*])([*])[*]')))
        if branch_points >= 3 and props['mw'] > 500:
            return True, "Highly branched macromolecular structure"

    # Special cases for known polymer types
    if any('mer' in comp.lower() for comp in smiles.split('.')):
        return True, "Contains explicit polymer notation"
        
    if num_components >= 3 and all(p['mw'] > 200 for p in component_props):
        return True, "Complex mixture of large molecules"

    return False, "Does not meet polymer characteristics"