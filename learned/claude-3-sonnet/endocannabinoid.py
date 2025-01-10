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
    Key classes include N-acylethanolamines and 2-acylglycerols.

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

    # Key structural patterns
    ethanolamine = Chem.MolFromSmarts("[NX3][CH2][CH2][OH]")
    glycerol_2_position = Chem.MolFromSmarts("[CH2X4][CHX4]([OH,OR])[CH2X4]")
    
    # Fatty acid chain patterns - more specific
    long_chain = Chem.MolFromSmarts("[CH2,CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    arachidonoyl = Chem.MolFromSmarts("CC=CC=CC=CC=CCC")  # Simplified pattern for multiple double bonds
    
    # Linkage patterns
    amide = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    ester = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    
    # Negative patterns - exclude these
    phosphate = Chem.MolFromSmarts("[$(P(=[O,S])([O,S])[O,S])]")
    multiple_esters = Chem.MolFromSmarts("([CX3](=[OX1])[OX2][#6])[CX3](=[OX1])[OX2][#6]")
    
    # Basic checks
    if mol.HasSubstructMatch(phosphate):
        return False, "Contains phosphate group (not characteristic of endocannabinoids)"
    
    if len(mol.GetSubstructMatches(ester)) > 2:
        return False, "Too many ester groups for endocannabinoid"
    
    # Count rings - endocannabinoids typically have 0-1 rings (epoxy possible)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 1:
        return False, "Too many rings for endocannabinoid"
        
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 30:
        return False, "Carbon count outside typical range for endocannabinoids"

    # Check for characteristic structural features
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine)
    has_glycerol = mol.HasSubstructMatch(glycerol_2_position)
    has_long_chain = mol.HasSubstructMatch(long_chain)
    has_arachidonoyl = mol.HasSubstructMatch(arachidonoyl)
    
    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    features = []
    
    # N-acylethanolamine type
    if has_ethanolamine:
        if not mol.HasSubstructMatch(amide):
            return False, "Ethanolamine present but missing required amide linkage"
        features.append("N-acylethanolamine")
        
    # 2-acylglycerol type
    elif has_glycerol:
        if not mol.HasSubstructMatch(ester):
            return False, "Glycerol present but missing required ester linkage"
        features.append("2-acylglycerol")
    else:
        return False, "Missing required ethanolamine or glycerol backbone"
    
    # Check fatty acid characteristics
    if not has_long_chain:
        return False, "Missing required long fatty acid chain"
    
    if has_arachidonoyl:
        features.append("arachidonoyl chain")
    
    if double_bonds > 0:
        features.append(f"{double_bonds} double bond{'s' if double_bonds > 1 else ''}")
    
    # Must have appropriate chain characteristics
    if not (has_arachidonoyl or (has_long_chain and double_bonds >= 1)):
        return False, "Fatty acid chain lacks characteristic unsaturation"
        
    # Additional structural requirements
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, "Molecular weight outside typical range for endocannabinoids"
        
    return True, f"Endocannabinoid: {', '.join(features)}"