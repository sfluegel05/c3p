"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:24913 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Most icosanoids should have around 20 carbons (allowing some variation due to metabolism)
    if c_count < 15 or c_count > 25:
        return False, f"Carbon count ({c_count}) outside typical range for icosanoids"
    
    # Should have multiple oxygens from oxidation
    if o_count < 2:
        return False, "Too few oxygen atoms for an icosanoid"

    # Look for carboxylic acid group (common in fatty acid derivatives)
    acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for carbon chain of suitable length
    chain_pattern = Chem.MolFromSmarts('[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]')
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No suitable carbon chain found"

    # Look for common icosanoid features
    features_found = []
    
    # Double bonds (common in fatty acid derivatives)
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    if mol.HasSubstructMatch(double_bond_pattern):
        features_found.append("double bonds")
    
    # Hydroxyl groups (from oxidation)
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    if mol.HasSubstructMatch(hydroxyl_pattern):
        features_found.append("hydroxyl groups")
    
    # Cyclopentane ring (as in prostaglandins)
    cyclopentane_pattern = Chem.MolFromSmarts('C1CCCC1')
    if mol.HasSubstructMatch(cyclopentane_pattern):
        features_found.append("cyclopentane ring")
        
    # Epoxide groups (as in EETs)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    if mol.HasSubstructMatch(epoxide_pattern):
        features_found.append("epoxide group")
    
    # Ketone groups (as in prostaglandins)
    ketone_pattern = Chem.MolFromSmarts('C(=O)C')
    if mol.HasSubstructMatch(ketone_pattern):
        features_found.append("ketone group")

    # Need at least some oxidation features to be an icosanoid
    if len(features_found) < 1:
        return False, "Lacks typical icosanoid oxidation features"

    # Calculate molecular weight - should be in reasonable range for C20 oxidized fatty acid
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for icosanoids"

    features_str = ", ".join(features_found)
    return True, f"Contains characteristic icosanoid features: {features_str}"