"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
A heterocyclic fatty acid containing an epoxide ring as part of its structure.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_epoxy_fatty_acid, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Check for epoxide ring
    epoxide = Chem.MolFromSmarts('[CR0]1[OR0][CR0]1')
    if not mol.HasSubstructMatch(epoxide):
        return False, "No epoxide ring found"
    
    # Count carbons and check chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon chain too short ({carbon_count} carbons) for fatty acid"
    if carbon_count > 30:
        return False, f"Carbon chain too long ({carbon_count} carbons) for fatty acid"
    
    # Validate fatty acid chain
    fatty_chain_patterns = [
        'CCCCCCCC',  # Minimum 8-carbon chain
        'CCCCC/C=C/CC',  # Unsaturated chain
        'CCCC/C=C/C/C=C/C',  # Polyunsaturated chain
        'CCCCC1OC1CC'  # Chain with epoxide
    ]
    
    chain_found = False
    for pattern in fatty_chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            chain_found = True
            break
    
    if not chain_found:
        return False, "No characteristic fatty acid chain found"
    
    # Count epoxide rings and validate their position
    epoxide_matches = len(mol.GetSubstructMatches(epoxide))
    if epoxide_matches > 3:
        return False, f"Too many epoxide rings ({epoxide_matches})"
        
    # Check epoxide is connected to main chain
    epoxide_chain = Chem.MolFromSmarts('C(CC)C1OC1C(CC)')
    if not mol.HasSubstructMatch(epoxide_chain):
        return False, "Epoxide not properly connected to main chain"
    
    # Count oxygen atoms (excluding carboxylic acid oxygens)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, f"Too few oxygen atoms ({oxygen_count})"
    if oxygen_count > 10:
        return False, f"Too many oxygen atoms ({oxygen_count})"
    
    # Check for non-fatty acid structures
    problematic_patterns = [
        '[#6]1[#6][#6]1',  # cyclopropane
        'C1CCCCC1',  # cyclohexane
        'C(=O)OC(=O)',  # anhydride
        'C(=O)NC(=O)',  # imide
        'C1CCOC1'  # tetrahydrofuran
    ]
    
    for pattern in problematic_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains non-fatty acid structural elements"
    
    # Count double bonds (excluding carboxylic acid)
    double_bond_count = rdMolDescriptors.CalcNumAliphicDoubleBonds(mol) - 1
    
    # Build description
    features = []
    if double_bond_count > 0:
        features.append(f"{double_bond_count} double bonds")
    features.append(f"{epoxide_matches} epoxide ring(s)")
    additional_oxygens = oxygen_count - (2 + epoxide_matches)  # Subtract COOH and epoxide oxygens
    if additional_oxygens > 0:
        features.append(f"{additional_oxygens} additional oxygen-containing groups")
    
    return True, f"Epoxy fatty acid with {carbon_count} carbons, " + ", ".join(features)