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

    # Basic chromanol/chroman core pattern
    # More permissive pattern that matches the core structure with various substitutions
    chromanol_pattern = """
        [#8]1-[#6](-[#6])(-[#6])-[#6]-[#6]-c2c1cc([#8])cc2
    """
    
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(chromanol_pattern.strip())):
        return False, "No chroman-6-ol core structure found"

    # Check for hydroxyl or modified hydroxyl (like esters) at position 6
    hydroxyl_patterns = [
        "[#8]1-[#6]-[#6]-[#6]-[#6]-c2c1cc([OH1])cc2",  # Free hydroxyl
        "[#8]1-[#6]-[#6]-[#6]-[#6]-c2c1cc([#8]-[#6])cc2"  # Modified hydroxyl
    ]
    
    hydroxyl_found = False
    for pattern in hydroxyl_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            hydroxyl_found = True
            break
            
    if not hydroxyl_found:
        return False, "Missing required hydroxyl group at position 6"

    # Check for the isoprenoid chain at position 2
    # More flexible patterns to catch both saturated and unsaturated chains
    chain_patterns = [
        # Saturated chain pattern
        "[#8]1-[#6](-[#6][#6][#6][#6])-[#6]-[#6]-c2c1cccc2",
        # Unsaturated chain pattern
        "[#8]1-[#6](-[#6][#6]=[#6][#6])-[#6]-[#6]-c2c1cccc2"
    ]
    
    chain_found = False
    for pattern in chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            chain_found = True
            break
            
    if not chain_found:
        return False, "Missing required chain at position 2"

    # Verify molecule size and composition
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 16:
        return False, "Insufficient carbons for tocol structure"
    
    if num_oxygens < 2:
        return False, "Insufficient oxygens for tocol structure"

    # Check for common modifications
    modifications = []
    
    # Check for ester modifications
    if mol.HasSubstructMatch(Chem.MolFromSmarts("O=C-O")):
        modifications.append("ester")
    
    # Check for carboxylic acid modifications
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O")):
        modifications.append("carboxylic acid")
    
    # Check for unsaturated bonds in chain
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        modifications.append("unsaturated chain")

    base_message = "Valid tocol structure with chroman-6-ol core and appropriate chain"
    if modifications:
        return True, f"{base_message} ({', '.join(modifications)})"
    
    return True, base_message