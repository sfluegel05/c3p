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

    # More comprehensive patterns for different flavonoid cores
    patterns = {
        # Basic chromane (flavan) core
        'chromane': '[#6]1[#6][#6][#6]2=[#6][#6]=[#6][#6]=[#6]2[O][#6]1',
        
        # Flavone/flavonol core with ketone
        'flavone': '[#6]1=[#6]-[#6](=[O])-[#6]2=[#6][#6]=[#6][#6]=[#6]2-[#8]-1',
        
        # Flavan-3-ol core (as in procyanidins)
        'flavan3ol': '[#6]1[#6][#6]([#8])[#6]2=[#6][#6]=[#6][#6]=[#6]2[#8][#6]1',
        
        # Common biflavonoid linkage patterns
        'c8_c8_link': 'c1c(c)c(O)c(c2)c1-c1c(c2)c(O)c',  # C8-C8 linkage
        'c8_c3_link': 'c1c(c)c(O)c(C-C2=COC)c1',  # C8-C3 linkage
        'c_o_c_link': 'c1cccc(Oc2cccc2)c1'  # C-O-C linkage
    }
    
    # Convert patterns to RDKit molecules
    pattern_mols = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}
    
    # Count core structures
    matches = {name: len(mol.GetSubstructMatches(pattern)) 
              for name, pattern in pattern_mols.items()}
    
    # Check for presence of basic requirements
    chromane_count = matches['chromane']
    flavone_count = matches['flavone']
    flavan3ol_count = matches['flavan3ol']
    
    total_cores = max(chromane_count, flavone_count + flavan3ol_count)
    
    # Check for linkage patterns
    has_c8_c8 = matches['c8_c8_link'] > 0
    has_c8_c3 = matches['c8_c3_link'] > 0
    has_c_o_c = matches['c_o_c_link'] > 0
    
    # Basic requirements checks
    if total_cores < 2:
        return False, f"Insufficient flavonoid cores (found {total_cores}, need at least 2)"
    
    # Count key features
    oh_pattern = Chem.MolFromSmarts('[OH]')
    oh_count = len(mol.GetSubstructMatches(oh_pattern))
    
    aromatic_rings = len(mol.GetAromaticRings())
    
    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450 or mol_wt > 900:
        return False, f"Molecular weight {mol_wt:.1f} outside typical biflavonoid range (450-900)"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Structural requirements for biflavonoids
    if aromatic_rings < 4:
        return False, f"Too few aromatic rings ({aromatic_rings}, need ≥4)"
    
    if oh_count < 2:
        return False, f"Too few hydroxyl groups ({oh_count}, need ≥2)"
    
    if c_count < 25 or c_count > 45:
        return False, f"Carbon count ({c_count}) outside typical range for biflavonoids"
    
    if o_count < 6:
        return False, f"Insufficient oxygen atoms ({o_count}, need ≥6)"
    
    # Check for linkage evidence
    if not (has_c8_c8 or has_c8_c3 or has_c_o_c):
        return False, "No characteristic biflavonoid linkage pattern found"
    
    # If all checks pass, classify as biflavonoid
    linkage_types = []
    if has_c8_c8: linkage_types.append("C8-C8")
    if has_c8_c3: linkage_types.append("C8-C3")
    if has_c_o_c: linkage_types.append("C-O-C")
    
    return True, f"Contains {total_cores} flavonoid cores with {', '.join(linkage_types)} linkage(s)"