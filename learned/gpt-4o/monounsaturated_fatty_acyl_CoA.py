"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated Fatty Acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for Coenzyme A (CoA) pattern - more flexible to handle variability
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A pattern found"

    # Find all C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=O.C=C")
    double_bond_matches = [match for match in mol.GetSubstructMatches(double_bond_pattern)]

    # Ensure that double bond counted is part of the fatty acyl chain separately
    if len(double_bond_matches) != 1:
        return False, f"Found {len(double_bond_matches)} carbon-carbon double bonds, need exactly 1 in the acyl chain"

    # Check for indicative acyl-CoA linkage
    acyl_linkage_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(acyl_linkage_pattern):
        return False, "No fatty acyl chain with proper linkage detected"
    
    return True, "Structure is a monounsaturated fatty acyl-CoA with one double bond in the acyl chain"