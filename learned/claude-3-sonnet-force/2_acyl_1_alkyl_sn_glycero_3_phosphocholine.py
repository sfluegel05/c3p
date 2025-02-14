"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:18035 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule belongs to the class 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
    based on its SMILES string. This class is defined as an alkyl,acyl-sn-glycero-3-phosphocholine
    with unspecified alkyl and acyl groups located at positions 1 and 2 respectively.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "No glycerol backbone found"

    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)([O-])OCC[N+](C)(C)C")
    phosphocholine_match = mol.GetSubstructMatch(phosphocholine_pattern)
    if not phosphocholine_match:
        return False, "No phosphocholine group found"

    # Check if phosphocholine group is attached to glycerol at position 3
    glycerol_atoms = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in glycerol_match]
    phosphocholine_atoms = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in phosphocholine_match]
    if glycerol_atoms[2] != "C" or phosphocholine_atoms[0] != "O":
        return False, "Phosphocholine group not attached to glycerol at position 3"

    # Check for alkyl group at position 1
    alkyl_atom = mol.GetAtomWithIdx(glycerol_match[0])
    alkyl_chain = []
    for neighbor in alkyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            alkyl_chain.append(neighbor.GetIdx())
            for atom in Chem.FindAllPathsOfLengthN(mol, alkyl_chain[-1], 8, useBonds=False):
                alkyl_chain.extend(atom)
    if len(alkyl_chain) < 4:
        return False, "Alkyl chain too short"

    # Check for acyl group at position 2
    acyl_atom = mol.GetAtomWithIdx(glycerol_match[1])
    acyl_neighbors = [neighbor for neighbor in acyl_atom.GetNeighbors() if neighbor.GetAtomicNum() == 8]
    if not acyl_neighbors:
        return False, "No acyl group found at position 2"
    acyl_oxygen = acyl_neighbors[0]
    acyl_carbon = [neighbor for neighbor in acyl_oxygen.GetNeighbors() if neighbor.GetAtomicNum() == 6][0]
    if acyl_carbon.GetTotalNumHs() != 0 or len(acyl_carbon.GetNeighbors()) != 3:
        return False, "Acyl group not found at position 2"

    return True, "Molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"