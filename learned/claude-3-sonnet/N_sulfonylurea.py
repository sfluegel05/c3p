"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    N-sulfonylureas contain a urea group where one nitrogen is connected to a sulfonyl group.
    General structure: R-SO2-NH-C(=O)-NH-R'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define exact N-sulfonylurea pattern:
    # [#6,#7,#8] - Carbon, Nitrogen or Oxygen (R group)
    # S(=O)(=O) - Sulfonyl group
    # [NH] - Nitrogen with at least one H
    # C(=O) - Carbonyl
    # [NH] - Second nitrogen with at least one H
    # [#6,#7,#8,#1] - R' group or H
    nsulf_pattern = Chem.MolFromSmarts('[#6,#7,#8]S(=O)(=O)[NH]C(=O)[NH][#6,#7,#8,#1]')
    
    # Alternative pattern for N-substituted variants
    nsulf_pattern2 = Chem.MolFromSmarts('[#6,#7,#8]S(=O)(=O)NC(=O)N[#6,#7,#8]')

    if not (mol.HasSubstructMatch(nsulf_pattern) or mol.HasSubstructMatch(nsulf_pattern2)):
        return False, "Missing required N-sulfonylurea substructure"

    # Exclude invalid cases
    
    # Pattern for sulfamate derivatives (O-SO2-NH2)
    sulfamate_pattern = Chem.MolFromSmarts('[O]S(=O)(=O)[NH2]')
    
    # Pattern for sulfonamide without urea (-SO2-NH-CO-)
    sulfonamide_pattern = Chem.MolFromSmarts('[#6,#7,#8]S(=O)(=O)[NH][CX3](=[O])[!N]')
    
    # Pattern for N-sulfonylcarbamate (R-O-CO-NH-SO2-R)
    sulfonylcarbamate_pattern = Chem.MolFromSmarts('[O][CH0]([!N])=O')

    if mol.HasSubstructMatch(sulfamate_pattern):
        return False, "Contains sulfamate group instead of N-sulfonylurea"
        
    if mol.HasSubstructMatch(sulfonylcarbamate_pattern):
        return False, "Contains N-sulfonylcarbamate instead of N-sulfonylurea"

    # Validate proper connectivity for each match
    for match in mol.GetSubstructMatches(nsulf_pattern):
        r1, s, n1, c, o, n2, r2 = match
        
        # Check oxidation states
        s_atom = mol.GetAtomWithIdx(s)
        if s_atom.GetTotalValence() != 6:
            continue  # Not a proper sulfonyl group
            
        # Verify carbonyl oxygen
        o_atom = mol.GetAtomWithIdx(o)
        if o_atom.GetTotalNumHs() != 0:
            continue  # Not a proper carbonyl
            
        # Verify urea nitrogens
        n1_atom = mol.GetAtomWithIdx(n1)
        n2_atom = mol.GetAtomWithIdx(n2)
        
        # At least one of the nitrogens should have a hydrogen
        if n1_atom.GetTotalNumHs() + n2_atom.GetTotalNumHs() == 0:
            continue
            
        return True, "Contains N-sulfonylurea group with proper connectivity"

    # Check alternative pattern
    for match in mol.GetSubstructMatches(nsulf_pattern2):
        r1, s, n1, c, o, n2, r2 = match
        
        s_atom = mol.GetAtomWithIdx(s)
        if s_atom.GetTotalValence() != 6:
            continue
            
        o_atom = mol.GetAtomWithIdx(o)
        if o_atom.GetTotalNumHs() != 0:
            continue
            
        return True, "Contains N-sulfonylurea group with proper connectivity"

    return False, "Structure lacks proper N-sulfonylurea connectivity"