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
        # Chromone core (more specific than benzopyran)
        'chromone': Chem.MolFromSmarts('O=C1CC(c2ccccc2)Oc2ccccc21'),
        
        # Flavan core
        'flavan': Chem.MolFromSmarts('C1CC(c2ccccc2)Oc2ccccc21'),
        
        # Common biflavonoid linkage patterns
        'c8_c8': Chem.MolFromSmarts('c1c(c2)c(O)cc(O)c1-c1c(c3)c(O)cc(O)c1'),
        'c3_c8': Chem.MolFromSmarts('c1c(c2)c(O)cc(O)c1-c1c(O)cc(O)c2c1'),
        'c6_c8': Chem.MolFromSmarts('c1c(c2)c(O)cc(-c3c(O)cc(O)c4c3)c1'),
        
        # Hydroxyl and methoxy groups
        'hydroxyl': Chem.MolFromSmarts('cO'),
        'methoxy': Chem.MolFromSmarts('cOC'),
        
        # Common oxygen bridge
        'o_bridge': Chem.MolFromSmarts('c1cccc(Oc2ccccc2)c1')
    }

    # Count rings
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = sum(1 for ring in ring_info.AtomRings() 
                             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_ring_count < 4:
        return False, f"Too few aromatic rings ({aromatic_ring_count}, need ≥4)"

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

    # Check for presence of hydroxyl/methoxy groups
    hydroxyl_count = len(mol.GetSubstructMatches(patterns['hydroxyl']))
    methoxy_count = len(mol.GetSubstructMatches(patterns['methoxy']))
    if hydroxyl_count + methoxy_count < 4:
        return False, f"Too few OH/OMe groups ({hydroxyl_count + methoxy_count}, need ≥4)"

    # Look for flavonoid cores
    chromone_matches = len(mol.GetSubstructMatches(patterns['chromone']))
    flavan_matches = len(mol.GetSubstructMatches(patterns['flavan']))
    
    # Check for common biflavonoid linkages
    c8_c8_links = len(mol.GetSubstructMatches(patterns['c8_c8']))
    c3_c8_links = len(mol.GetSubstructMatches(patterns['c3_c8']))
    c6_c8_links = len(mol.GetSubstructMatches(patterns['c6_c8']))
    o_bridges = len(mol.GetSubstructMatches(patterns['o_bridge']))
    
    total_cores = chromone_matches + flavan_matches
    if total_cores < 1:
        return False, "No flavonoid cores found"
        
    total_linkages = c8_c8_links + c3_c8_links + c6_c8_links + o_bridges
    if total_linkages == 0:
        return False, "No characteristic biflavonoid linkage found"

    # Additional structural checks
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 6:
        return False, f"Too few rings ({ring_count}, need ≥6)"

    # If all checks pass, it's likely a biflavonoid
    return True, f"Contains flavonoid cores with characteristic biflavonoid linkage(s) and {hydroxyl_count} OH groups"