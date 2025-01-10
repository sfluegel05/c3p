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

    # Check for classic icosanoid patterns (prostaglandins, thromboxanes)
    cyclopentane = Chem.MolFromSmarts("C1CCCC1")
    if mol.HasSubstructMatch(cyclopentane):
        return False, "Contains cyclopentane ring (characteristic of classic icosanoids)"

    # Count oxygenated functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]")  # hydroxyl group
    epoxy_pattern = Chem.MolFromSmarts("[OX2r3]")    # epoxy oxygen
    
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    epoxy_matches = mol.GetSubstructMatches(epoxy_pattern)
    
    hydroxy_count = len(hydroxy_matches) - 1  # Subtract 1 for carboxylic acid OH
    epoxy_count = len(epoxy_matches)
    
    if hydroxy_count == 0 and epoxy_count == 0:
        return False, "No hydroxyl or epoxy groups found"

    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    
    if double_bond_count == 0:
        return False, "No double bonds found"

    # Check for conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    conjugated_count = len(mol.GetSubstructMatches(conjugated_pattern))

    # Count total oxygens (should have multiple for oxygenation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for oxygenated fatty acid"

    # Build feature description
    features = []
    if hydroxy_count > 0:
        features.append(f"{hydroxy_count} hydroxyl groups")
    if epoxy_count > 0:
        features.append(f"{epoxy_count} epoxy groups")
    if conjugated_count > 0:
        features.append(f"conjugated double bond system")
    else:
        features.append(f"{double_bond_count} double bonds")

    # Validate molecular weight (should be around 300-400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} outside expected range"

    return True, "Nonclassic icosanoid: " + ", ".join(features)