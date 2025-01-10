"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible substructure pattern for backbone
    glycerol_pattern = Chem.MolFromSmarts("OCCO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Define phosphoserine pattern, avoiding strict stereochemistry
    phosphoserine_pattern = Chem.MolFromSmarts("COP(=O)(O)OC[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Look for acyl group O=C-C-O connectivity at the sn-1 position
    acyl_sn1_pattern = Chem.MolFromSmarts("O-C(=O)-")
    if not mol.HasSubstructMatch(acyl_sn1_pattern):
        return False, "Missing acyl group at sn-1 position"

    # Calculate Length of the Acyl Chain to affirm typical length
    acyl_carbon_chain_pattern = Chem.MolFromSmarts("C(=O)OC[C;R0]") # Flexible C chain from ester
    matches = mol.GetSubstructMatches(acyl_carbon_chain_pattern)
    
    found_long_chain = any(
        rdMolDescriptors.CalcNumRotatableBonds(Chem.PathToSubmol(mol, match)) >= 10
        for match in matches
    )
    
    if not found_long_chain:
        return False, "Acyl chain too short for 1-acyl-sn-glycero-3-phosphoserine"

    return True, "Molecule contains 1-acyl-sn-glycero-3-phosphoserine"