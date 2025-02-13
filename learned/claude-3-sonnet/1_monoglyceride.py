"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at position 1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_1_monoglyceride, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Look for exactly two free hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) != 2:
        return False, f"Found {len(hydroxyl_matches)} free hydroxyl groups, need exactly 2"

    # Pattern for 1-monoglyceride specific arrangement
    # [CH2] with ester, connected to [CH] with OH, connected to [CH2] with OH
    mono_1_pattern = Chem.MolFromSmarts("[CH2X4][OX2][CX3](=[OX1])-[!$([OH])]")
    if not mol.HasSubstructMatch(mono_1_pattern):
        return False, "Ester group not at position 1"

    # Verify carbon chain length (should have at least 4 carbons total including glycerol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short"

    # Count oxygens (should have exactly 4: 2 from ester, 2 from hydroxyls)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 4:
        return False, "Must have exactly 4 oxygens"

    # Check for reasonable molecular weight (should be > 100 Da for smallest 1-monoglyceride)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for 1-monoglyceride"

    return True, "Contains glycerol backbone with one ester bond at position 1 and two free hydroxyl groups"