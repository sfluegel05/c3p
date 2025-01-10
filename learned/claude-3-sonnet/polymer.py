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
from collections import Counter

def has_regular_repeating_units(mol, min_repeats=3):
    """Helper function to detect regular repeating structural units"""
    # Common polymer repeat unit patterns
    repeat_patterns = [
        # Prenyl/isoprene units
        ('CC(C)=CCC', 'isoprene'),
        # Common polymer linkages
        ('COC', 'ether'),
        ('CC(=O)O', 'ester'),
        ('CN(C)C(=O)', 'amide'),
        ('CSC', 'thioether'),
        ('CNC(=O)', 'peptide'),
        # Sugar units
        ('C1C(O)C(O)C(O)C(O)C1', 'sugar'),
    ]
    
    for pattern, unit_type in repeat_patterns:
        matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
        if matches >= min_repeats:
            return True, unit_type
    return False, None

def analyze_components(components):
    """Helper function to analyze molecular components"""
    mols = []
    for comp in components:
        mol = Chem.MolFromSmiles(comp)
        if mol is None:
            continue
        mols.append({
            'mol': mol,
            'mw': Descriptors.ExactMolWt(mol),
            'charge': Chem.GetFormalCharge(mol),
            'atoms': mol.GetNumAtoms(),
            'rings': rdMolDescriptors.CalcNumRings(mol)
        })
    return mols

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    components = smiles.split('.')
    if not components:
        return False, "Invalid SMILES string"
        
    # Analyze components
    comp_data = analyze_components(components)
    if not comp_data:
        return False, "Invalid SMILES string"
    
    # Check for prenyl polymers (special case)
    for comp in comp_data:
        mol = comp['mol']
        isoprene_pattern = Chem.MolFromSmarts('CC(C)=CCC')
        if isoprene_pattern and len(mol.GetSubstructMatches(isoprene_pattern)) >= 4:
            return True, "Contains multiple prenyl/isoprene repeat units"
    
    # Check for regular repeating units in each component
    for comp in comp_data:
        has_repeats, unit_type = has_regular_repeating_units(comp['mol'])
        if has_repeats:
            return True, f"Contains regular repeating {unit_type} units"
    
    # Analyze component mixture characteristics
    num_components = len(comp_data)
    if num_components >= 2:
        # Check for charged species (salt forms)
        charges = [c['charge'] for c in comp_data]
        identical_components = len(set(components)) < len(components)
        
        if any(c != 0 for c in charges) and identical_components:
            return True, "Ionic polymer with repeating components"
            
        # Check for identical large components
        large_comps = [c for c in comp_data if c['mw'] > 400]
        if len(large_comps) >= 2 and identical_components:
            return True, "Contains multiple identical large components"
    
    # Check for very large molecules with internal repetition
    for comp in comp_data:
        mol = comp['mol']
        if comp['mw'] > 1000:
            # Look for internal symmetry and repetition
            fragments = Chem.GetMolFrags(mol, asMols=True)
            if len(fragments) > 1:
                fragment_smiles = [Chem.MolToSmiles(f) for f in fragments]
                if len(set(fragment_smiles)) < len(fragment_smiles):
                    return True, "Large molecule with internal repeating units"
    
    return False, "Does not meet polymer characteristics"