"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine has a glycerol backbone with acyl substituents at the 1- and
    2-hydroxy positions and a phosphoserine group at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key features:
    # 1. Glycerol backbone with phosphoserine and ester links at 1 and 2 positions
    glycerol_phosphoserine_pattern = Chem.MolFromSmarts(
        "[C@@H](CO[P]([O-])([O])=O)O[C@H](CO[C@H](N)C(=O)O)OC(=O)[C]"
    )
    
    if not mol.HasSubstructMatch(glycerol_phosphoserine_pattern):
        return False, "No matching glycerol-phosphoserine backbone found, or incorrect acyl chain positions"
    
    # 2. Pattern for acyl chains connected via ester bonds
    ester_chain_pattern = Chem.MolFromSmarts(
        "C(=O)O[C@@H](COP([O-])([O])=O)OC(=O)C"
    )

    if not mol.HasSubstructMatch(ester_chain_pattern):
        return False, "No ester bonds recognized between glycerol backbone and acyl chains"

    # Overall match confirmation - if any substructure needed was missing, function would have returned by now
    return True, "Structure matches criteria for 3-sn-phosphatidyl-L-serine"