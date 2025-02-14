"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (C-C-C with 2 esters and 1 phosphate)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2][CX3](=[OX1]))[CH2X4][OX2][P]")
    if not mol.HasSubstructMatch(glycerol_pattern):
         return False, "No glycerol-phosphate backbone found"

    # Check for 2 ester groups (fatty acid chains)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"


    # Check for a myo-inositol ring (six membered ring with 5 OH groups)
    inositol_pattern = Chem.MolFromSmarts("C1[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring found"

    # Check for at least one phosphate group attached to the inositol ring
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C]([O])[OP]")
    inositol_phosphate_matches = mol.GetSubstructMatches(inositol_phosphate_pattern)
    if len(inositol_phosphate_matches) < 1:
        return False, "No phosphate group attached to the inositol ring found"

    return True, "Contains glycerol backbone, two fatty acid chains, a phosphate group linked to glycerol, a myo-inositol ring, and at least one phosphate group on the inositol ring"