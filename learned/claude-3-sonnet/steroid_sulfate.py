"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # More flexible steroid core patterns to match different variations
    steroid_patterns = [
        # Basic steroid core with flexible bond types
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1",
        # Alternative pattern allowing for aromatic rings
        "c1cc2[C,c]~[C,c]~[C,c]3~[C,c]~[C,c]~[C,c]4~[C,c]~[C,c]~[C,c]~[C,c]~4~[C,c]~[C,c]3~[C,c]~[C,c]2c1",
        # Pattern for estrone-like steroids
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6,c]~[#6,c]~[#6,c]~[#6,c]~[#6,c]~3~[#6]~[#6]~2~[#6]~1"
    ]
    
    found_steroid = False
    for pattern in steroid_patterns:
        steroid_pat = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(steroid_pat):
            found_steroid = True
            break
            
    if not found_steroid:
        return False, "No steroid core structure found"
    
    # Various sulfate group patterns
    sulfate_patterns = [
        "OS(=O)(=O)O",  # neutral sulfate
        "OS(=O)(=O)[O-]",  # deprotonated sulfate
        "[#6]-[#8]S(=O)(=O)[O-]",  # alternative form
        "[#6]-[#8]S(=O)(=O)O",  # alternative neutral form
        "[#6]-OS([O-])(=O)=O"  # another ionic form
    ]
    
    found_sulfate = False
    for pattern in sulfate_patterns:
        sulfate_pat = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(sulfate_pat):
            found_sulfate = True
            break
            
    if not found_sulfate:
        return False, "No sulfate group found"
    
    # Count sulfate groups
    total_sulfates = sum(len(mol.GetSubstructMatches(Chem.MolFromSmarts(pat))) 
                        for pat in sulfate_patterns)
    
    # Verify sulfate is attached to carbon (part of steroid)
    sulfate_linkage_patterns = [
        "[#6]-[#8]-[#16](=[#8])(=[#8])-[#8]",
        "[#6]-[#8]S([O-])(=O)=O"
    ]
    
    proper_linkage = False
    for pattern in sulfate_linkage_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            proper_linkage = True
            break
            
    if not proper_linkage:
        return False, "Sulfate group not properly connected to steroid core"
    
    # Basic size check for steroid structure
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 17:  # Most steroids have 17+ carbons
        return False, "Too few carbons for steroid structure"
    
    if total_sulfates == 1:
        return True, "Found steroid core with one sulfate group"
    else:
        return True, f"Found steroid core with {total_sulfates} sulfate groups"