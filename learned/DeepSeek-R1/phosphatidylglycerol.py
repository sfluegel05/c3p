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
    A phosphatidylglycerol consists of a glycerol backbone with two fatty acid esters
    and a phosphatidyl group (glycerol connected via phosphate).

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

    # Check for exactly two ester groups (O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Core structure: Phosphate connecting two glycerol units
    # Check for two glycerol backbones connected via phosphate
    # Pattern for glycerol (three consecutive carbons with oxygen attachments)
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-O)-[CH2]")
    phosphate_pattern = Chem.MolFromSmarts("[PX4](-O)(-O)(=O)-O")

    # Find all phosphate atoms connected to two glycerol units
    found = False
    for p in p_atoms:
        # Get neighboring oxygens connected via single bonds
        o_neighbors = [n for n in p.GetNeighbors() if n.GetAtomicNum() == 8 and 
                      mol.GetBondBetweenAtoms(p.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
        if len(o_neighbors) < 2:
            continue

        # Check each oxygen is part of a glycerol backbone
        glycerol_count = 0
        for o in o_neighbors:
            # Get the carbon attached to this oxygen
            c_neighbors = [n for n in o.GetNeighbors() if n.GetAtomicNum() == 6]
            if not c_neighbors:
                continue
            c = c_neighbors[0]
            # Check if this carbon is part of a glycerol structure
            if mol.GetSubstructMatch(glycerol_pattern):
                glycerol_count += 1

        if glycerol_count >= 2:
            found = True
            break

    if not found:
        return False, "Glycerol-phosphate-glycerol backbone not found"

    # Verify one glycerol has two esters (head) and the other has hydroxyl (tail)
    # Check ester groups are attached to glycerol backbones
    # Get all ester groups
    esters = [mol.GetAtomWithIdx(a) for match in ester_matches for a in match if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]
    # Check esters are attached to glycerol carbons
    glycerol_carbons = set()
    for match in mol.GetSubstructMatches(glycerol_pattern):
        glycerol_carbons.update(match)
    ester_attached = False
    for ester in esters:
        for neighbor in ester.GetNeighbors():
            if neighbor.GetIdx() in glycerol_carbons:
                ester_attached = True
                break
    if not ester_attached:
        return False, "Esters not attached to glycerol backbone"

    # Check for at least one hydroxyl in the tail glycerol
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups in tail glycerol"

    # Check molecular weight (typical PG > 600)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    # Check rotatable bonds for fatty acid chains
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds < 10:
        return False, f"Insufficient rotatable bonds ({rot_bonds}) for fatty acids"

    return True, "Contains glycerol-phosphate-glycerol backbone with two fatty acid esters"