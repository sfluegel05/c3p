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
    
    # More flexible steroid core pattern that can match different oxidation states
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Check for carboxylic acid group or its conjugates
    acid_patterns = [
        "[CX3](=O)[OX2H1]",  # carboxylic acid
        "[CX3](=O)NCC(=O)[OH]",  # glycine conjugate
        "[CX3](=O)NCCS(=O)(=O)[OH]",  # taurine conjugate
        "[CX3](=O)[O-]",  # carboxylate
    ]
    
    has_acid = False
    for pattern in acid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_acid = True
            break
            
    if not has_acid:
        return False, "No carboxylic acid group or conjugates found"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count ({c_count}) outside typical range for bile acids"
    
    if o_count < 3:
        return False, "Insufficient oxygen atoms for bile acid"
    
    # Look for hydroxyl or oxo groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    oxo_pattern = Chem.MolFromSmarts("[#6]=O")
    
    hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    oxo_groups = len(mol.GetSubstructMatches(oxo_pattern))
    
    if hydroxyls + oxo_groups < 1:
        return False, "No hydroxyl or oxo groups found"
    
    # Molecular features check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 700:
        return False, f"Molecular weight ({mol_wt}) outside typical range for bile acids"
    
    # Ring analysis
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"
    
    # Check for sufficient sp3 carbons (characteristic of steroid core)
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_carbons < 12:
        return False, "Too few sp3 carbons for bile acid structure"
        
    # Look for characteristic side chain pattern
    side_chain = Chem.MolFromSmarts("[CH2][CH2][CH2]C(=O)[OH]")
    if not mol.HasSubstructMatch(side_chain):
        side_chain2 = Chem.MolFromSmarts("[CH2][CH2]C(=O)[OH]")
        if not mol.HasSubstructMatch(side_chain2):
            return False, "Missing characteristic bile acid side chain"

    return True, "Matches bile acid structure with steroid core, appropriate functional groups, and characteristic features"