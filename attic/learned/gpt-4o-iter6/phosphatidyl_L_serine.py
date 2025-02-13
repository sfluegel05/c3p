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

    # Phosphate group pattern: P=O with two oxygens, one being part of a glycerol-like structure
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
        
    # Serine esterified pattern - amino acid backbone (serine) pattern
    serine_pattern = Chem.MolFromSmarts("C[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine group esterified to the phosphate"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone pattern found"

    # Check for two ester groups connected to the glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups; found {len(ester_matches)}, need at least 2"

    return True, "Molecule meets the criteria for phosphatidyl-L-serine"

# Examples can be tested using this function with the SMILES strings provided.