"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: CHEBI:17515 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol is a glycerol derivative where one primary hydroxy group
    is replaced by a phosphatidyl group (glycerol connected via phosphate).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of phosphorus (phosphate group)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphate group found"

    # Define core structure pattern: glycerol connected via phosphate to another glycerol
    # Pattern: glycerol (C-O-P-O-C) where each C is part of a glycerol backbone
    core_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-O-[P](=O)(-O-)-O-)-[CH2]."  # First glycerol part
        "O-[P](=O)(-O-)-O-[CH2]-[CH](-O)-[CH2]"   # Phosphate connecting to second glycerol
    )
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Glycerol-phosphate-glycerol backbone not found"

    # Check for at least two ester groups (O-C=O) attached to the first glycerol
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester groups, need at least 2"

    # Verify esters are attached to the first glycerol backbone
    first_glycerol = Chem.MolFromSmarts("[CH2]-[CH](-O-[P])-[CH2]")
    first_glyc_matches = mol.GetSubstructMatches(first_glycerol)
    if not first_glyc_matches:
        return False, "First glycerol backbone not found"
    
    # Collect carbons in first glycerol
    glyc_carbons = set()
    for match in first_glyc_matches:
        glyc_carbons.update(match[:3])  # First three atoms are C-C-C in glycerol

    # Check if at least two esters are attached to the first glycerol
    ester_count = 0
    for ester in ester_matches:
        oxygen_idx = ester[0]
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        # Check adjacent carbon is part of the first glycerol
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in glyc_carbons:
                ester_count +=1
                break

    if ester_count < 2:
        return False, f"Only {ester_count} ester(s) attached to first glycerol"

    # Check the second glycerol has at least one hydroxyl group
    second_glycerol = Chem.MolFromSmarts("[CH2]-[CH](-O)-[CH2]")
    second_glyc_matches = mol.GetSubstructMatches(second_glycerol)
    if not second_glyc_matches:
        return False, "Second glycerol backbone not found"
    
    hydroxyl_found = False
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    for match in second_glyc_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0:
                hydroxyl_found = True
                break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "Second glycerol lacks hydroxyl group"

    return True, "Contains glycerol-phosphate-glycerol backbone with at least two fatty acid esters"