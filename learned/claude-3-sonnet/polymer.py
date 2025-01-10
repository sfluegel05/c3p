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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains multiple components (separated by dots)
    components = smiles.count('.')
    
    # Calculate molecular descriptors
    mol_wt = Descriptors.ExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Look for repeating units using common patterns
    repeating_carbon_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2]")
    repeating_matches = len(mol.GetSubstructMatches(repeating_carbon_chain))
    
    # Check for presence of metal ions
    metal_pattern = Chem.MolFromSmarts("[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,In,Sn,Pb,Bi]")
    has_metals = mol.HasSubstructMatch(metal_pattern) if metal_pattern else False
    
    # Define criteria for polymer classification
    is_large = mol_wt > 500
    has_many_atoms = num_atoms > 30
    has_long_chains = num_rotatable_bonds > 10
    has_repeating_units = repeating_matches >= 2
    has_multiple_components = components >= 1
    has_multiple_rings = num_rings >= 2

    # Classification logic
    if is_large and (has_repeating_units or has_long_chains):
        return True, "Contains repeating units or long chains with high molecular weight"
        
    if has_multiple_components and has_many_atoms:
        return True, "Multiple component mixture with large molecular size"
        
    if has_metals and has_multiple_components:
        return True, "Metal-containing polymer complex"
        
    if has_repeating_units and has_multiple_rings and has_long_chains:
        return True, "Complex structure with repeating units and multiple rings"
        
    if mol_wt > 1000 and (has_long_chains or has_multiple_rings):
        return True, "Very large molecule with complex structure"
        
    if components >= 2 and has_many_atoms and (has_multiple_rings or has_long_chains):
        return True, "Multi-component mixture with complex molecular structure"

    # Count carbon chain length
    carbon_chain = Chem.MolFromSmarts("C~C~C~C~C~C")
    if carbon_chain and len(mol.GetSubstructMatches(carbon_chain)) > 2:
        return True, "Contains multiple long carbon chains"

    # If none of the above criteria are met, it's probably not a polymer
    return False, "Does not meet polymer characteristics"