"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: fatty alcohols - aliphatic alcohols with carbon chains of 3+ atoms
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alcohol groups (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    
    if not alcohol_matches:
        return False, "No aliphatic alcohol group found"
    
    # Count carbons in longest chain
    # First, get all carbons that are not part of aromatic rings
    aliphatic_carbons = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            aliphatic_carbons.add(atom.GetIdx())
    
    if len(aliphatic_carbons) < 3:
        return False, "Carbon chain too short (minimum 3 carbons required)"
        
    # Check if OH is part of carboxylic acid
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "OH group is part of carboxylic acid"
    
    # Check for ester pattern
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]")
    if mol.HasSubstructMatch(ester_pattern):
        # If there's an ester, make sure there's still a free alcohol group
        # that's not part of the ester
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        ester_oxygens = set(match[3] for match in ester_matches)
        alcohol_oxygens = set(match[1] for match in alcohol_matches)
        
        if not (alcohol_oxygens - ester_oxygens):
            return False, "Contains ester but no free alcohol group"
    
    # Additional checks for reasonable molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:  # Arbitrary upper limit for reasonable fatty alcohols
        return False, "Molecular weight too high for typical fatty alcohol"
        
    # Count total carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count > 40:  # Arbitrary upper limit
        return False, "Carbon count too high for typical fatty alcohol"
        
    # Success case - construct detailed reason
    features = []
    features.append(f"Contains {len(alcohol_matches)} alcohol group(s)")
    features.append(f"Has {c_count} total carbons")
    
    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    triple_bonds = len(mol.GetSubstructMatches(triple_bond_pattern))
    
    if double_bonds > 0 or triple_bonds > 0:
        features.append(f"Unsaturated ({double_bonds} double bonds, {triple_bonds} triple bonds)")
    else:
        features.append("Saturated")
    
    # Check for branching
    branched = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
            branched = True
            break
    features.append("Branched" if branched else "Unbranched")
    
    return True, "; ".join(features)