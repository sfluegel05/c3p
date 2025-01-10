"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are signaling lipids that activate cannabinoid receptors.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for key structural features
    # Ethanolamine group (-NCCO)
    ethanolamine_pattern = Chem.MolFromSmarts("[NX3][CH2][CH2][OH]")
    # Glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    # Ester linkage
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    # Amide linkage
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    # Ether linkage
    ether_pattern = Chem.MolFromSmarts("[#6][OX2][#6]")
    # Double bonds
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    # Long aliphatic chain (at least 4 carbons)
    aliphatic_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2]")
    
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    has_long_chain = mol.HasSubstructMatch(aliphatic_chain)

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:  # Endocannabinoids typically have 0-1 rings
        return False, "Too many rings for endocannabinoid"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 30:
        return False, "Carbon count outside typical range for endocannabinoids"

    # Must have long aliphatic chain
    if not has_long_chain:
        return False, "Missing required long aliphatic chain"

    # Must have either ethanolamine or glycerol backbone
    if not (has_ethanolamine or has_glycerol):
        return False, "Missing required ethanolamine or glycerol group"
    
    # Must have either ester, amide, or ether linkage
    if not (has_ester or has_amide or has_ether):
        return False, "Missing required ester, amide, or ether linkage"
        
    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:  # Tightened molecular weight range
        return False, "Molecular weight outside typical range for endocannabinoids"

    # Build feature list
    features = []
    if has_ethanolamine:
        features.append("ethanolamine")
    if has_glycerol:
        features.append("glycerol")
    if has_ester:
        features.append("ester")
    if has_amide:
        features.append("amide")
    if has_ether:
        features.append("ether")
    if double_bonds >= 1:
        features.append(f"{double_bonds} double bond{'s' if double_bonds > 1 else ''}")
    
    # Classification logic:
    # Must have:
    # 1. Long aliphatic chain
    # 2. Either ethanolamine or glycerol
    # 3. Either ester, amide, or ether linkage
    # 4. Appropriate molecular weight
    # 5. Not too many rings
    
    return True, f"Endocannabinoid containing: {', '.join(features)}"