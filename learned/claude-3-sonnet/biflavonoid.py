"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid consists of two flavonoid units joined by a single bond or atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavonoid core patterns
    patterns = {
        # Benzopyran core (more general)
        'benzopyran': Chem.MolFromSmarts('c1ccc2OC(C)CCc2c1'),
        
        # Flavone core
        'flavone': Chem.MolFromSmarts('O=C1CC(c2ccccc2)Oc2ccccc21'),
        
        # Flavan core
        'flavan': Chem.MolFromSmarts('C1CC(c2ccccc2)Oc2ccccc21'),
        
        # Hydroxyl group
        'hydroxyl': Chem.MolFromSmarts('cO'),
        
        # Common biaryl linkage
        'biaryl': Chem.MolFromSmarts('c-c'),
        
        # Ether linkage
        'ether': Chem.MolFromSmarts('cOc')
    }

    # Count aromatic rings
    aromatic_rings = len(mol.GetAromaticRings())
    if aromatic_rings < 4:
        return False, f"Too few aromatic rings ({aromatic_rings}, need ≥4)"

    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450 or mol_wt > 900:
        return False, f"Molecular weight {mol_wt:.1f} outside typical biflavonoid range (450-900)"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:
        return False, f"Too few carbons ({c_count}, need ≥25)"
    if o_count < 6:
        return False, f"Too few oxygens ({o_count}, need ≥6)"

    # Check for presence of basic structural features
    hydroxyl_count = len(mol.GetSubstructMatches(patterns['hydroxyl']))
    if hydroxyl_count < 2:
        return False, f"Too few hydroxyl groups ({hydroxyl_count}, need ≥2)"

    # Look for flavonoid cores
    benzopyran_matches = len(mol.GetSubstructMatches(patterns['benzopyran']))
    flavone_matches = len(mol.GetSubstructMatches(patterns['flavone']))
    flavan_matches = len(mol.GetSubstructMatches(patterns['flavan']))
    
    total_cores = max(benzopyran_matches, flavone_matches + flavan_matches)
    if total_cores < 2:
        return False, f"Insufficient flavonoid cores (found {total_cores}, need ≥2)"

    # Check for linkage between units
    biaryl_links = len(mol.GetSubstructMatches(patterns['biaryl']))
    ether_links = len(mol.GetSubstructMatches(patterns['ether']))
    
    if biaryl_links == 0 and ether_links == 0:
        return False, "No characteristic biflavonoid linkage found"

    # Additional structural checks
    ring_systems = rdMolDescriptors.CalcNumRings(mol)
    if ring_systems < 6:
        return False, f"Too few ring systems ({ring_systems}, need ≥6)"

    # If all checks pass, it's likely a biflavonoid
    linkage_desc = []
    if biaryl_links > 0:
        linkage_desc.append("direct C-C")
    if ether_links > 0:
        linkage_desc.append("ether")
    
    return True, f"Contains {total_cores} flavonoid cores connected by {' and '.join(linkage_desc)} linkage(s)"