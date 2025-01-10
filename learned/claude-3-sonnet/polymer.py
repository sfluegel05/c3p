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

def analyze_ionic_complex(components):
    """Helper function to analyze ionic complexes"""
    mols = []
    total_charge = 0
    has_counterions = False
    
    for comp in components:
        mol = Chem.MolFromSmiles(comp)
        if mol is None:
            continue
        charge = Chem.GetFormalCharge(mol)
        total_charge += charge
        if abs(charge) > 0:
            has_counterions = True
        mols.append({
            'mol': mol,
            'charge': charge,
            'mw': Descriptors.ExactMolWt(mol)
        })
    
    # Check for charge balance and presence of counterions
    return has_counterions and abs(total_charge) <= 1

def has_polymer_characteristics(mol):
    """Helper function to detect polymer-specific characteristics"""
    # Check molecular weight
    mw = Descriptors.ExactMolWt(mol)
    if mw < 400:  # Too small to be a polymer unit
        return False
        
    # Count rotatable bonds (flexibility characteristic of polymers)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 8:  # Too rigid for typical polymer
        return False
        
    # Check for branching
    branching = rdMolDescriptors.CalcNumBranches(mol)
    if branching < 2:  # Linear small molecules
        return False
        
    return True

def has_repeating_units(mol, min_repeats=4):
    """Helper function to detect meaningful repeating structural units"""
    # Polymer-specific repeat unit patterns
    polymer_patterns = [
        # Prenyl/isoprene units (with spacing)
        ('CC(C)=CCC[R]CC(C)=C', 'prenyl'),
        # Common polymer backbones
        ('CCCCCCC', 'alkyl chain'),
        ('COC(=O)CCCC', 'polyester'),
        ('CN(C)C(=O)CCC', 'polyamide'),
        ('CSC(C)CC', 'polysulfide'),
        ('CC(O)C(O)CC', 'polyol'),
    ]
    
    for pattern, unit_type in polymer_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol and len(mol.GetSubstructMatches(pattern_mol)) >= min_repeats:
            return True, unit_type
    return False, None

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
    
    # Check for ionic polymers/complexes
    if len(components) > 1 and analyze_ionic_complex(components):
        return True, "Ionic polymer/complex with balanced charges"
    
    # Analyze each component
    for comp in components:
        mol = Chem.MolFromSmiles(comp)
        if mol is None:
            continue
            
        # Check for polymer characteristics
        if has_polymer_characteristics(mol):
            # Look for repeating units
            has_repeats, unit_type = has_repeating_units(mol)
            if has_repeats:
                return True, f"Contains characteristic polymer structure with {unit_type} repeating units"
                
        # Special case: Check for large prenyl/terpene polymers
        isoprene_pattern = Chem.MolFromSmarts('CC(C)=CCC')
        if isoprene_pattern and len(mol.GetSubstructMatches(isoprene_pattern)) >= 5:
            return True, "Natural polymer with multiple isoprene units"
            
    # Check for multi-component systems
    if len(components) >= 3:
        # Look for repeating components
        comp_counter = Counter(components)
        if any(count >= 2 for count in comp_counter.values()):
            return True, "Multi-component polymer system with repeating units"
            
    return False, "Does not meet polymer characteristics"