"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:78512 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are biologically active signalling molecules made by 
    oxygenation of C20 fatty acids, excluding classic icosanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Must have exactly 20 carbons, found {c_count}"

    # Count oxygens (should have at least 3 - one from COOH plus at least two more)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for nonclassic icosanoid"

    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    conjugated_double_bond_pattern = Chem.MolFromSmarts("C=CC=C")
    
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    conjugated_matches = len(mol.GetSubstructMatches(conjugated_double_bond_pattern))
    
    if double_bond_matches < 2:
        return False, "Insufficient double bonds"
    
    # Look for common oxygen-containing functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OH1]")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern))
    epoxy_count = len(mol.GetSubstructMatches(epoxy_pattern))
    
    if hydroxy_count == 0 and epoxy_count == 0:
        return False, "No hydroxyl or epoxy groups found"

    # Look for patterns that would indicate classic icosanoids
    prostaglandin_pattern = Chem.MolFromSmarts("[CH2][CH2][CH]1[CH2][CH2][CH]([CH2][CH2][CH2]C(=O)[OH])[CH]1")
    if mol.HasSubstructMatch(prostaglandin_pattern):
        return False, "Contains prostaglandin core structure"

    # Check for long carbon chain
    carbon_chain = Chem.MolFromSmarts("CCCCCC")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chain found"

    # Calculate molecular weight (should be around 300-400)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 450:
        return False, "Molecular weight outside expected range"

    # If all checks pass, it's likely a nonclassic icosanoid
    features = []
    if conjugated_matches > 0:
        features.append(f"{conjugated_matches} conjugated double bond systems")
    if hydroxy_count > 0:
        features.append(f"{hydroxy_count} hydroxyl groups")
    if epoxy_count > 0:
        features.append(f"{epoxy_count} epoxy groups")
    
    return True, f"C20 oxygenated fatty acid with {', '.join(features)}"