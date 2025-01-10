"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate is a glycerol backbone with a phosphate group at position 3,
    and a single acyl group attached at either position 1 or position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with phosphate at position 3
    glycerol3p_pattern = Chem.MolFromSmarts("O[P](=O)(O)[O][C@H](CO)CO")
    if not mol.HasSubstructMatch(glycerol3p_pattern):
        return False, "No glycerol 3-phosphate backbone found"

    # Check for ester-linked acyl group at position 1 or 2
    # Define ester linkage pattern attached to glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No acyl ester group found"
    elif len(ester_matches) > 1:
        return False, "More than one acyl ester group found"

    # Check that the acyl group is attached to position 1 or 2
    # Search for acyl group attached to glycerol backbone
    acyl_positions = []
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[2])  # The carbon atom connected to the ester oxygen
        neighbors = [nbr.GetIdx() for nbr in ester_atom.GetNeighbors()]
        # Look for connection to glycerol backbone carbons
        for idx in neighbors:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and atom.GetDegree() == 3:
                acyl_positions.append(idx)

    if len(acyl_positions) == 0:
        return False, "Acyl group not attached to glycerol backbone"
    elif len(acyl_positions) > 1:
        return False, "Multiple acyl groups attached to glycerol backbone"

    # Ensure phosphate group is at position 3
    # Define phosphate group pattern attached to glycerol
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)(O)[O][C@H](CO)CO")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not attached at position 3"

    # Check for additional acyl groups
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) > 1:
        return False, "Additional acyl groups found"

    # Confirm stereochemistry (sn configuration)
    # For simplicity, we can check if chiral centers are specified
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if len(chiral_centers) == 0:
        return False, "No chiral centers found, stereochemistry not specified"

    return True, "Contains glycerol 3-phosphate backbone with a single acyl group attached at position 1 or 2"