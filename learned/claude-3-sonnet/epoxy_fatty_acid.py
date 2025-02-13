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
from rdkit.Chem.AllChem import FindMolChiralCenters

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
    epoxide = Chem.MolFromSmarts('[C]1[O][C]1')
    if not mol.HasSubstructMatch(epoxide):
        return False, "No epoxide ring found"
    
    # Count carbons and check chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon chain too short ({carbon_count} carbons) for fatty acid"
    if carbon_count > 30:
        return False, f"Carbon chain too long ({carbon_count} carbons) for fatty acid"
    
    # Count epoxide rings
    epoxide_matches = len(mol.GetSubstructMatches(epoxide))
    if epoxide_matches > 3:
        return False, f"Too many epoxide rings ({epoxide_matches})"
    
    # Count oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, f"Too few oxygen atoms ({oxygen_count})"
    if oxygen_count > 10:
        return False, f"Too many oxygen atoms ({oxygen_count})"
    
    # Check for excessive rings (excluding epoxide)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count > epoxide_matches + 1:  # Allow epoxide rings plus maybe one other
        return False, "Contains too many ring systems"
    
    # Count double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    alkene_count = double_bond_count - 1  # Subtract carboxylic C=O
    
    # Check for fatty acid chain patterns - multiple possibilities
    fatty_patterns = [
        'CCCCCC',  # straight chain
        'CCC(C)CC',  # branched
        'CC=CC',  # unsaturated
        'CCCCC1OC1',  # epoxide containing
    ]
    
    chain_found = False
    for pattern in fatty_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            chain_found = True
            break
            
    if not chain_found:
        return False, "No characteristic fatty acid chain found"
    
    # Check for characteristic epoxy fatty acid patterns
    epoxy_patterns = [
        'CC1OC1CC(=O)O',  # terminal epoxide
        'CC1OC1CCC',  # mid-chain epoxide
        'C=CC1OC1CC',  # unsaturated epoxide
    ]
    
    characteristic_pattern = False
    for pattern in epoxy_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            characteristic_pattern = True
            break
    
    if not characteristic_pattern:
        return False, "No characteristic epoxy fatty acid pattern found"
    
    # Build description
    features = []
    if alkene_count > 0:
        features.append(f"{alkene_count} double bonds")
    features.append(f"{epoxide_matches} epoxide ring(s)")
    if oxygen_count > 3:
        features.append(f"{oxygen_count-3} additional oxygen-containing groups")
    
    return True, f"Epoxy fatty acid with {carbon_count} carbons, " + ", ".join(features)