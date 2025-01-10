"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids occurring in bile, typically present
    as sodium salts of their amides with glycine or taurine.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific steroid core patterns that can match different oxidation states
    # Pattern includes both saturated and unsaturated versions of the core
    steroid_patterns = [
        # Basic steroid core with flexible bond types
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6]4[#6][#6]3[#6]2[#6]1",
        # Core with potential ketone groups
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6](=O)[#6]4[#6][#6]3[#6]2[#6]1",
        # Core with potential double bonds
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6]=,:[#6][#6]4[#6][#6][#6][#6]4[#6][#6]3[#6]2[#6]1"
    ]
    
    has_core = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
            
    if not has_core:
        return False, "No suitable steroid core found"

    # Check for carboxylic acid group or its conjugates
    acid_patterns = [
        "[CX3](=O)[OX2H1]",  # carboxylic acid
        "[CX3](=O)[O-]",  # carboxylate
        "[CX3](=O)NCC(=O)[OH]",  # glycine conjugate
        "[CX3](=O)NCCS(=O)(=O)[OH]"  # taurine conjugate
    ]
    
    has_acid = False
    for pattern in acid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_acid = True
            break
            
    if not has_acid:
        return False, "No carboxylic acid group or conjugates found"

    # Common hydroxyl/oxo positions in bile acids
    hydroxyl_patterns = [
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6]4[#6][#6]3[#6]2[#6]1[OX2H1]",  # 3-OH
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][OX2H1][#6]4[#6][#6][#6][#6]4[#6][#6]3[#6]2[#6]1",  # 7-OH
        "[#6]1[#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][OX2H1][#6][#6]4[#6][#6]3[#6]2[#6]1"   # 12-OH
    ]
    
    hydroxyl_count = 0
    for pattern in hydroxyl_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            hydroxyl_count += 1

    # Count total hydroxyls and oxo groups
    total_hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    total_oxo_pattern = Chem.MolFromSmarts("[#6]=O")
    
    total_hydroxyls = len(mol.GetSubstructMatches(total_hydroxyl_pattern))
    total_oxo = len(mol.GetSubstructMatches(total_oxo_pattern))
    
    if total_hydroxyls + total_oxo < 1:
        return False, "Insufficient hydroxyl/oxo groups for bile acid"

    # Basic molecular properties
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (20 <= c_count <= 30):
        return False, f"Carbon count ({c_count}) outside typical range for bile acids"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (300 <= mol_wt <= 700):
        return False, f"Molecular weight ({mol_wt}) outside typical range for bile acids"

    # Ring count check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Side chain check - allow for variations
    side_chain_patterns = [
        "[CH2][CH2][CH2]C(=O)[OH]",  # standard
        "[CH2][CH2]C(=O)[OH]",       # shortened
        "[CH2][CH2][CH2]C(=O)[O-]",  # ionized
        "[CH2][CH2][CH2]C(=O)N"      # conjugated
    ]
    
    has_side_chain = False
    for pattern in side_chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_side_chain = True
            break
            
    if not has_side_chain:
        return False, "Missing characteristic bile acid side chain"

    return True, "Matches bile acid structure with appropriate steroid core, functional groups, and characteristic features"