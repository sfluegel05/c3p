"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: CHEBI:17855 diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride has a glycerol backbone with two fatty acid chains via ester bonds,
    and the third hydroxyl group is either free or substituted with an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (three carbons in a row with appropriate bonds)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check exactly two ester groups attached to the glycerol backbone
    ester_attached_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[CX3]=[OX1]")
    ester_attached_matches = mol.GetSubstructMatches(ester_attached_pattern)
    if len(ester_attached_matches) != 2:
        return False, f"Found {len(ester_attached_matches)} ester groups attached to glycerol, need exactly 2"

    # Find the third oxygen in the glycerol (not part of ester groups)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    backbone_indices = glycerol_matches[0]  # Take first match
    backbone_carbons = [mol.GetAtomWithIdx(idx) for idx in backbone_indices]

    # Collect all oxygens attached to the backbone carbons
    glycerol_oxygens = []
    for carbon in backbone_carbons:
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                glycerol_oxygens.append(neighbor)

    # Get oxygen indices from ester groups (match[1] is the ester oxygen in our SMARTS pattern)
    ester_oxygen_indices = {match[1] for match in ester_attached_matches}
    
    # Identify the third oxygen not part of the ester groups
    third_oxygens = [o for o in glycerol_oxygens if o.GetIdx() not in ester_oxygen_indices]
    if len(third_oxygens) != 1:
        return False, "Could not identify third oxygen in glycerol"
    third_oxygen = third_oxygens[0]

    # Check third oxygen is either -OH or O-alkyl
    if third_oxygen.GetDegree() == 1:
        pass  # Valid -OH
    elif third_oxygen.GetDegree() == 2:
        # Check if connected to a carbon (alkyl group)
        other_neighbors = [n for n in third_oxygen.GetNeighbors() if n not in backbone_carbons]
        if not other_neighbors or other_neighbors[0].GetAtomicNum() != 6:
            return False, "Third oxygen is O-alkyl but not connected to carbon"
    else:
        return False, "Third oxygen has invalid bonding"

    # Check oxygen count (3 from glycerol + 2 from ester carbonyls)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 5:
        return False, f"Expected 5 oxygens, found {o_count}"

    # Check molecular weight (typical diglycerides >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} < 300)"

    # Check carbon count (at least 15)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons ({c_count})"

    # Check rotatable bonds (indicates chain length)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, f"Too few rotatable bonds ({n_rotatable})"

    return True, "Contains glycerol backbone with two ester-linked fatty acids and one hydroxyl/O-alkyl group"