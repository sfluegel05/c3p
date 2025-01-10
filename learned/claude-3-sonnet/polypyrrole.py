"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a polypyrrole and reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define pyrrole pattern: five-membered ring with one nitrogen
    # [nH] represents pyrrole nitrogen, #1 assigns it number 1
    # The rest defines carbons in the ring with single/double bonds
    pyrrole_pattern = Chem.MolFromSmarts('[nH,n;R1]1:[c,C;R1]:[c,C;R1]:[c,C;R1]:[c,C;R1]:1')
    
    # Find all pyrrole matches
    matches = mol.GetSubstructMatches(pyrrole_pattern)
    num_pyrroles = len(matches)
    
    if num_pyrroles == 0:
        return False, "No pyrrole rings found"
    elif num_pyrroles == 1:
        return False, "Only one pyrrole ring found, need at least two"
    
    # Additional check for porphyrin-like structures (four connected pyrroles)
    porphyrin_pattern = Chem.MolFromSmarts('[nH,n]1:[c,C]:[c,C]:[c,C]:[c,C]:1-[c,C]1:[c,C]:[c,C]2:[nH,n]:[c,C]:[c,C]:[c,C]:[c,C]:2-[c,C]2:[c,C]:[c,C]3:[nH,n]:[c,C]:[c,C]:[c,C]:[c,C]:3-[c,C]3:[c,C]:[c,C]4:[nH,n]:[c,C]:[c,C]:[c,C]:[c,C]:4-[c,C]:3:[c,C]:2:[c,C]:1')
    
    is_porphyrin = mol.HasSubstructMatch(porphyrin_pattern)
    
    # Check for metal coordination
    metal_pattern = Chem.MolFromSmarts('[Fe,Mg,Zn,Co,Pd,B]')
    has_metal = mol.HasSubstructMatch(metal_pattern)
    
    # Construct detailed reason
    reason = f"Contains {num_pyrroles} pyrrole rings"
    if is_porphyrin:
        reason += " arranged in a porphyrin-like structure"
    if has_metal:
        reason += " with metal coordination"
        
    return True, reason