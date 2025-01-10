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

    # More flexible steroid core patterns
    steroid_patterns = [
        # Basic cyclopentanoperhydrophenanthrene core with flexible connectivity
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~2~[#6]~1",
        # Alternative pattern allowing for ketones and modified rings
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6,#8]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~2~[#6]~1",
        # Pattern for systems with double bonds
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]=,:[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~2~[#6]~1",
        # More generic pattern for modified systems
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6,#8]~[#6,#8]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~2~[#6]~1"
    ]
    
    has_core = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_core = True
            break
            
    if not has_core:
        return False, "No suitable steroid core found"

    # Expanded carboxylic acid patterns
    acid_patterns = [
        "[CX3](=O)[OX2H1,OX1-]",  # acid or carboxylate
        "[CX3](=O)[NX3][CH2][CX3](=O)[OX2H1,OX1-]",  # glycine conjugate
        "[CX3](=O)[NX3][CH2][CH2][SX4](=O)(=O)[OX2H1,OX1-]",  # taurine conjugate
        "[CX3](=O)[OX2][#6]",  # ester (might be prodrug or derivative)
    ]
    
    has_acid = False
    for pattern in acid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_acid = True
            break
            
    if not has_acid:
        return False, "No carboxylic acid group or conjugates found"

    # Expanded hydroxyl/oxo patterns for bile acids
    functional_group_patterns = [
        # Common hydroxyl positions (3,7,12)
        "[#6][#6]([#6])[OX2H1]",  # general hydroxyl
        "[#6][CX3](=O)[#6]",      # ketone
        "[#6][#6]([#6])=[O]"      # aldehyde
    ]
    
    functional_group_count = 0
    for pattern in functional_group_patterns:
        matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
        functional_group_count += matches

    if functional_group_count < 1:
        return False, "Insufficient functional groups for bile acid"

    # Basic molecular properties check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (20 <= c_count <= 35):  # Expanded range to include conjugates
        return False, f"Carbon count ({c_count}) outside typical range for bile acids"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (300 <= mol_wt <= 800):  # Expanded range to include conjugates
        return False, f"Molecular weight ({mol_wt}) outside typical range for bile acids"

    # Ring count and connectivity check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Side chain patterns - more flexible
    side_chain_patterns = [
        "[#6][#6][#6]C(=O)[O,N]",     # standard chain
        "[#6][#6]C(=O)[O,N]",         # shortened chain
        "[#6][#6][#6][#6]C(=O)[O,N]", # extended chain
        "[#6][#6][#6]C(=O)",          # partial match
    ]
    
    has_side_chain = False
    for pattern in side_chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_side_chain = True
            break
            
    if not has_side_chain:
        return False, "Missing characteristic bile acid side chain"

    return True, "Matches bile acid structure with appropriate steroid core, functional groups, and characteristic features"