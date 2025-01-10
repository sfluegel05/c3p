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
    # Double bonds
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    # Epoxy group
    epoxy_pattern = Chem.MolFromSmarts("[#6][OX2][#6]")

    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Basic requirements
    if not (has_ethanolamine or has_glycerol):
        return False, "Missing required ethanolamine or glycerol group"
    
    if not (has_ester or has_amide):
        return False, "Missing required ester or amide linkage"
    
    # Must have sufficient carbon chain length
    if c_count < 15:
        return False, "Carbon chain too short for endocannabinoid"
        
    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Most endocannabinoids are lipids with MW between 300-400 Da
    if mol_wt < 250 or mol_wt > 800:
        return False, "Molecular weight outside typical range for endocannabinoids"

    # Classify based on structural features
    features = []
    if has_ethanolamine:
        features.append("ethanolamine")
    if has_glycerol:
        features.append("glycerol")
    if has_ester:
        features.append("ester")
    if has_amide:
        features.append("amide")
    if double_bonds >= 1:
        features.append(f"{double_bonds} double bonds")
    if has_epoxy:
        features.append("epoxy group")
    
    # Return classification
    if len(features) >= 2:  # Must have at least 2 characteristic features
        return True, f"Endocannabinoid containing: {', '.join(features)}"
    else:
        return False, "Insufficient characteristic features for endocannabinoid classification"