"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: CHEBI:67194 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are characterized by long hydrocarbon chains connected to
    ethanolamide, glycerol, or similar groups via ester, amide, or ether bonds.

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

    # Check for ethanolamide group (-NCCO)
    ethanolamide_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][OX2H]")
    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)

    # Check for glycerol backbone (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)

    # Check for ether linkage (C-O-C)
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    has_ether = mol.HasSubstructMatch(ether_pattern)

    if not (has_ethanolamide or has_glycerol or has_ether):
        return False, "No ethanolamide, glycerol, or ether group found"

    # Check for long hydrocarbon chain (at least 10 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check for ester, amide, or ether linkage
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern) or has_ether):
        return False, "No ester, amide, or ether linkage found"

    # Check molecular weight (typically >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for endocannabinoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 16:
        return False, "Too few carbons for endocannabinoid"
    if o_count < 2:
        return False, "Too few oxygens for endocannabinoid"

    return True, "Contains long hydrocarbon chain connected to ethanolamide, glycerol, or similar group via ester, amide, or ether bond"