"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone with one acyl group (esterified fatty acid)
    and the remaining two positions occupied by either hydrogen or alkyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible glycerol backbone pattern (C-C-C with at least 1 oxygen)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check if ester is attached to glycerol backbone
    ester_attached = False
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in ester_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                if mol.HasSubstructMatch(glycerol_pattern):
                    ester_attached = True
                    break
    if not ester_attached:
        return False, "Ester group not attached to glycerol backbone"

    # Check for acyl group (minimum 4 carbons)
    acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 1:
        return False, "Missing acyl group"

    # Check molecular weight - monoacylglycerols typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for monoacylglycerol"

    # Count hydroxyl groups (should be at least 1)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0)
    if hydroxyl_count < 1:
        return False, "Not enough hydroxyl groups"

    return True, "Contains glycerol backbone with one acyl group and two H/alkyl groups"