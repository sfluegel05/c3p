"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine consists of a glycerol backbone with two fatty acid esters,
    a phosphate group, and an L-serine moiety connected via the phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (3 carbons with at least 2 oxygen connections)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for two ester groups (fatty acids)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for L-serine moiety with correct stereochemistry
    # Pattern: O-C-C@H(N)-C(=O)O (allowing for protonation states)
    serine_pattern = Chem.MolFromSmarts("[OX2]-[CH2]-[C@H]([NH2,NH3+])-[CX3](=O)[OX2H1,OX1-]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "Missing L-serine moiety with correct stereochemistry"

    # Verify serine's oxygen is connected to a phosphate group
    has_serine_phosphate = False
    for match in mol.GetSubstructMatches(serine_pattern):
        # Oxygen is first atom in the SMARTS pattern [O]-C-...
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Check if connected to phosphorus
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'P':
                has_serine_phosphate = True
                break
        if has_serine_phosphate:
            break
    
    if not has_serine_phosphate:
        return False, "Serine moiety not connected to phosphate group"

    # Check for phosphate group presence (P with >=3 oxygen neighbors)
    phosphorus = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus:
        return False, "No phosphorus atom found"
    
    valid_phosphate = False
    for p_atom in phosphorus:
        oxygen_count = sum(1 for neighbor in p_atom.GetNeighbors() if neighbor.GetAtomicNum() == 8)
        if oxygen_count >= 3:
            valid_phosphate = True
            break
    if not valid_phosphate:
        return False, "Invalid phosphate group configuration"

    # Check fatty acid chain length (at least 6 carbons each)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, "Insufficient long-chain substituents for fatty acids"

    return True, "Contains glycerol backbone with two fatty acid esters and L-serine-connected phosphate"