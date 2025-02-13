"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is composed of two phosphatidic acid molecules linked to a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count key atoms
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Must have exactly 2 phosphorus atoms
    if phosphorus_count != 2:
        return False, f"Must have exactly 2 phosphorus atoms, found {phosphorus_count}"

    # Look for phosphate groups (-P(=O)(O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O,OH])[O,OH]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"
    
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches != 2:
        return False, f"Found {phosphate_matches} phosphate groups, need exactly 2"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches != 4:
        return False, f"Found {ester_matches} ester groups, need exactly 4 for fatty acid chains"

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = len(mol.GetSubstructMatches(glycerol_pattern))
    if glycerol_matches < 2:
        return False, "Missing required glycerol backbones"

    # Check for fatty acid chains (long carbon chains)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = len(mol.GetSubstructMatches(fatty_acid_pattern))
    if fatty_acid_matches < 4:
        return False, "Missing required fatty acid chains"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 20:
        return False, "Too few rotatable bonds for cardiolipin"

    # Check molecular weight - cardiolipins typically >1200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for cardiolipin"

    # Count carbons - cardiolipins typically have >60 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 60:
        return False, f"Too few carbons ({carbon_count}) for cardiolipin"

    # Verify connectivity between phosphates and glycerol
    phosphoglycerol_pattern = Chem.MolFromSmarts("[P](=[O])([O,OH])[OCH2][CH][CH2]")
    if len(mol.GetSubstructMatches(phosphoglycerol_pattern)) < 2:
        return False, "Missing required phosphoglycerol connectivity"

    return True, "Contains two phosphatidic acids linked to glycerol with four fatty acid chains"