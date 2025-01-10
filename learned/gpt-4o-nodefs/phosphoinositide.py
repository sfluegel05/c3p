"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide typically contains a myo-inositol ring with one or more phosphate groups
    and fatty acid chains.

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
    
    # Define SMARTS pattern for myo-inositol ring
    inositol_pattern = Chem.MolFromSmarts("C1(OC)C(O)C(O)C(O)C(O)O1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol structure found"

    # Check for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "Less than one phosphate group found"
    
    # Optional: Check for connected phosphates (common in derivatives)
    # connected_phosphate_pattern = Chem.MolFromSmarts("[O,P](=O)(O)[O;R0]")
    
    # Check for fatty acid tails via ester linkage to glycerol backbones
    # Look for ester linkage as part of glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CO)CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester-linked glycerol backbone found"

    # Reinterpretation of long hydrocarbon chain through more flexible approach
    # to accommodate various chain lengths seen in examples
    long_chain_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCC)C") # Ideal: At minimum, detect one chain.
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)

    if len(long_chain_matches) < 1:
        return False, "No sufficient long hydrocarbon chains detected"
    
    return True, "Compound matches a phosphoinositide with typical features"

# Testing example structure
# smiles_str = "O=C(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)O[C@@H](COP(O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)OP(O)(=O)O)O)O)(=O)O)COC(=O)CCCCCCCCCCCCCCCCC"
# print(is_phosphoinositide(smiles_str))