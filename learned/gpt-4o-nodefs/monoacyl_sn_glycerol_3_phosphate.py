"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify glycerol structure (C3 chain with OH groups similar to C[OH]C[OH]C[OH])
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Identify the presence of a single acyl chain (long carbon chain with ester linkage to glycerol)
    acyl_pattern = Chem.MolFromSmarts("C(=O)OCC")  # Simplified ester linkage pattern
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Expected 1 acyl chain, found {len(acyl_matches)}"
    
    # Ensure a phosphate group is present at the third position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")  # Simplified phosphate group attached to glycerol
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found attached to glycerol"

    return True, "Contains glycerol, single acyl chain, and phosphate group at the third position"