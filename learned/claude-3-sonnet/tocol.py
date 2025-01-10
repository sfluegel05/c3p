"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: tocol compounds
Definition: A chromanol with a chroman-6-ol skeleton substituted at position 2 
by a saturated or triply-unsaturated hydrocarbon chain of three isoprenoid units
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_tocol, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic chromanol core pattern - more flexible than before
    # Matches the basic 6-hydroxychroman skeleton with possible methyl substitutions
    chromanol_core = """
        [#6]1-[#6]-[#6]-c2c(C1)c([OH0,OH1,O&$(OC(=O))])c([#6,H])c([#6,H])c2
    """
    
    core_pattern = Chem.MolFromSmarts(chromanol_core.strip())
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No chromanol core structure found"

    # Check for oxygen heterocycle in the correct position
    chroman_oxygen = """
        [#8]1-[#6](-[#6,#1])(-[#6,#1])-[#6]-[#6]-c2c1c([#6,H])c([OH0,OH1,O&$(OC(=O))])c([#6,H])c2
    """
    oxygen_pattern = Chem.MolFromSmarts(chroman_oxygen.strip())
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "Incorrect oxygen position in chromanol core"

    # Check for proper substitution at position 2
    # More flexible pattern that allows for various chain types
    position2_pattern = """
        [#8]1-[#6](-[#6][#6,#1])-[#6]-[#6]-c2c1c([#6,H])c([OH0,OH1,O&$(OC(=O))])c([#6,H])c2
    """
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(position2_pattern.strip())):
        return False, "Incorrect substitution at position 2"

    # Check for isoprenoid chain - both saturated and unsaturated versions
    # More flexible patterns that capture various modifications
    isoprenoid_patterns = [
        # Saturated chain
        "[#6]-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]",
        # Unsaturated chain
        "[#6]-[#6](-[#6])=[#6]-[#6]-[#6](-[#6])=[#6]-[#6]-[#6](-[#6])=[#6]",
        # Modified chain (e.g., with terminal COOH)
        "[#6]-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]-[#6]-[#6](-[#6])-[#6]-[#6](-[O,C])",
    ]

    chain_found = False
    for pattern in isoprenoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            chain_found = True
            break

    if not chain_found:
        return False, "Missing required isoprenoid chain"

    # Check for common modifications that are allowed
    modifications = []
    
    # Check for ester modifications
    ester_pattern = Chem.MolFromSmarts("[OH0,OH1,O&$(OC(=O))]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) > 1:
        modifications.append("ester derivative")

    # Check for carboxylic acid modifications
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if mol.HasSubstructMatch(acid_pattern):
        modifications.append("carboxylic acid")

    # Basic requirements check
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if carbon_count < 16:
        return False, "Insufficient carbon count for tocol structure"
    
    if oxygen_count < 2:
        return False, "Insufficient oxygen count for tocol structure"

    base_message = "Valid tocol structure with chromanol core and appropriate chain"
    if modifications:
        return True, f"{base_message} ({', '.join(modifications)})"
    
    return True, base_message