"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Classifies a molecule as a phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Correcting the phosphate pattern to match phosphate ester found in phosphatidyl-L-serine
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(OCC[N])(OC)OC(=O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphatidyl linkage found"
        
    # Correcting the serine linkage: Capture the serine esterified to phosphate and glycerol
    serine_pattern = Chem.MolFromSmarts("N[C@@H](C(=O)O)COP(=O)(O)OC")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No esterified serine group pattern found"

    # Glycerol backbone bound through ester linkages to two fatty acids
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H](COP(=O)(O)O)OC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester linkages found"

    # Verify two ester groups bound (indicative of fatty acid chains)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups; found {len(ester_matches)}, need at least 2"

    return True, "Molecule meets the criteria for phosphatidyl-L-serine"