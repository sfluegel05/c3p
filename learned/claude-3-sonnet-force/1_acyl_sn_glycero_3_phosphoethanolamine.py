"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule belongs to the 1-acyl-sn-glycero-3-phosphoethanolamine class based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the (R)-configuration of the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](OP(O)([O-,+0])OCC[NH3+])([CH2][OH])[CH2][OH]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No (R)-glycerophosphate backbone found"

    # Check for the presence of at least one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Check for the presence of carbon chains (not necessarily fatty acids)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "No carbon chains found"

    # Check for the presence of an aminoethyl group
    aminoethyl_pattern = Chem.MolFromSmarts("[NH2][CX4][CX4][OX2]")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl group found"

    return True, "Molecule belongs to the 1-acyl-sn-glycero-3-phosphoethanolamine class"