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
    # Three-membered ring containing oxygen
    epoxide = Chem.MolFromSmarts('[C]1[O][C]1')
    if not mol.HasSubstructMatch(epoxide):
        return False, "No epoxide ring found"
    
    # Count carbons to verify fatty acid chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:  # Most fatty acids have 12+ carbons
        return False, f"Carbon chain too short ({carbon_count} carbons) for fatty acid"
        
    # Verify the molecule is not too large for a fatty acid
    if carbon_count > 30:  # Most fatty acids have less than 30 carbons
        return False, f"Carbon chain too long ({carbon_count} carbons) for fatty acid"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:  # Typical range for epoxy fatty acids
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for epoxy fatty acids"
    
    # Count epoxide rings
    epoxide_matches = len(mol.GetSubstructMatches(epoxide))
    if epoxide_matches > 3:  # Most epoxy fatty acids have 1-3 epoxide rings
        return False, f"Too many epoxide rings ({epoxide_matches})"
        
    # Check for reasonable number of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:  # Minimum: 1 for epoxide + 2 for carboxylic acid
        return False, f"Too few oxygen atoms ({oxygen_count})"
    if oxygen_count > 8:  # Maximum: typically not more than 8
        return False, f"Too many oxygen atoms ({oxygen_count})"
    
    # Calculate degree of unsaturation
    double_bonds = rdMolDescriptors.CalcNumDoubleBonds(mol)
    if double_bonds > 6:  # Most epoxy fatty acids have 0-6 double bonds
        return False, f"Too many double bonds ({double_bonds})"
        
    # Verify the presence of an aliphatic chain
    aliphatic_chain = Chem.MolFromSmarts('[CH2][CH2][CH2]')
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No aliphatic chain found"
    
    return True, f"Contains epoxide ring and fatty acid chain with {carbon_count} carbons"