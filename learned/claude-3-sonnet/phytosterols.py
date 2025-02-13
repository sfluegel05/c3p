"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and vary 
    only in carbon side chains and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - more flexible than before
    # Four fused rings with some flexibility in bond types
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1")
    if steroid_core is None:
        return False, "Error in steroid core SMARTS pattern"
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for hydroxyl group
    hydroxyl = Chem.MolFromSmarts("[OX2H1]")
    if hydroxyl is None:
        return False, "Error in hydroxyl SMARTS pattern"
        
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "No hydroxyl group found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Phytosterols typically have 27-30 carbons, but we'll be a bit flexible
    if c_count < 25 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical range for phytosterols (25-35)"

    # Most phytosterols have 1-3 oxygen atoms (more if glycosylated)
    if o_count < 1 or o_count > 8:
        return False, f"Oxygen count ({o_count}) outside typical range for phytosterols"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 800:  # Increased range to account for glycosides
        return False, "Molecular weight outside typical range for phytosterols"

    # Look for branched aliphatic side chain
    side_chain = Chem.MolFromSmarts("[CH2,CH3][CH,C]([CH3,CH2])[CH2,CH][CH2,CH]")
    if side_chain is None:
        return False, "Error in side chain SMARTS pattern"
        
    if not mol.HasSubstructMatch(side_chain):
        return False, "No characteristic phytosterol side chain found"

    # Check for characteristics that would exclude it from being a phytosterol
    allowed_elements = {1, 6, 8} # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Contains elements other than C, H, and O"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Calculate number of double bonds
    double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bonds > 5:
        return False, "Too many double bonds for a phytosterol"

    return True, "Contains steroid core with hydroxyl group and characteristic phytosterol structure"