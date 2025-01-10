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
    if phosphorus_count != 2:
        return False, f"Must have exactly 2 phosphorus atoms, found {phosphorus_count}"

    # Look for the complete cardiolipin backbone structure:
    # Two phosphatidic acids (each with glycerol + phosphate) connected to a central glycerol
    cardiolipin_backbone = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]OP(=O)(O)OC[CH]CO")
    if not mol.HasSubstructMatch(cardiolipin_backbone):
        return False, "Missing required cardiolipin backbone structure"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches != 4:
        return False, f"Found {ester_matches} ester groups, need exactly 4 for fatty acid chains"

    # Look for fatty acid chains (long carbon chains)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = len(mol.GetSubstructMatches(fatty_acid_pattern))
    if fatty_acid_matches < 4:
        return False, "Missing required long fatty acid chains"

    # Check molecular weight - cardiolipins typically >1200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for cardiolipin"

    # Count carbons - cardiolipins typically have >60 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 60:
        return False, f"Too few carbons ({carbon_count}) for cardiolipin"

    # Count oxygens - cardiolipins should have at least 13 oxygens 
    # (4 esters + 2 phosphates + 3 glycerol)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 13:
        return False, f"Too few oxygens ({oxygen_count}) for cardiolipin"

    # Verify the complete connectivity pattern:
    # Central glycerol with two phosphate groups each connected to another glycerol
    complete_pattern = Chem.MolFromSmarts(
        "[CH2X4][CHX4][CH2X4]OP(=O)(O)OC[CH]([CH2]OC(=O)[#6])[CH2]OC(=O)[#6]"
    )
    if not mol.HasSubstructMatch(complete_pattern):
        return False, "Missing required cardiolipin connectivity pattern"

    return True, "Contains two phosphatidic acids linked to glycerol with four fatty acid chains"