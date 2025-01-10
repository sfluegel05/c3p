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

    # Key structural patterns
    ethanolamine_pattern = Chem.MolFromSmarts("[NX3][CH2][CH2][OH]")
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    fatty_acid_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]")  # At least 6 carbons
    conjugated_double_bonds = Chem.MolFromSmarts("[CH2]~[CH1]=[CH1]~[CH2]~[CH1]=[CH1]")
    
    # Linkage patterns
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    ether_pattern = Chem.MolFromSmarts("[#6][OX2][#6]")
    
    # Check for key features
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_fatty_chain = mol.HasSubstructMatch(fatty_acid_chain)
    has_conjugated_db = mol.HasSubstructMatch(conjugated_double_bonds)
    
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    
    # Count rings - endocannabinoids can have epoxy rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 3:  # Allow up to 3 rings (some have epoxy rings)
        return False, "Too many rings for endocannabinoid"

    # Must have either ethanolamine or glycerol backbone
    if not (has_ethanolamine or has_glycerol):
        return False, "Missing required ethanolamine or glycerol group"
    
    # Must have fatty acid chain
    if not has_fatty_chain:
        return False, "Missing required fatty acid chain"
    
    # Must have proper linkage
    if not (has_ester or has_amide or has_ether):
        return False, "Missing required ester, amide, or ether linkage"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Build feature list for classification reason
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
    if double_bonds > 0:
        features.append(f"{double_bonds} double bond{'s' if double_bonds > 1 else ''}")
    if has_conjugated_db:
        features.append("conjugated double bonds")
    
    # Additional structural requirements
    if has_ethanolamine:
        # For ethanolamine derivatives, require amide linkage
        if not has_amide:
            return False, "Ethanolamine derivatives must have amide linkage"
    
    if has_glycerol:
        # For glycerol derivatives, require ester or ether linkage
        if not (has_ester or has_ether):
            return False, "Glycerol derivatives must have ester or ether linkage"
    
    # At least one of these conditions must be true for typical endocannabinoids
    if not any([
        double_bonds >= 3,  # Multiple double bonds like in arachidonoyl derivatives
        (has_fatty_chain and (has_ethanolamine or has_glycerol)),  # Long chain with proper head group
        has_conjugated_db  # Conjugated double bond system
    ]):
        return False, "Structure lacks characteristic endocannabinoid features"

    return True, f"Endocannabinoid containing: {', '.join(features)}"