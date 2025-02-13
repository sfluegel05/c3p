"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: Monounsaturated Fatty Acyl-CoA
"""
from rdkit import Chem

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

    # Look for Coenzyme A (CoA) pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A pattern found"

    # Identify the acyl chain which starts from the thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S[C;X4]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found, implying absence of acyl chain"

    # Check for exactly one C=C double bond in the acyl chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    # Filter double bonds to ensure they're part of the fatty acyl chain
    acyl_chain_bonds = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0 and 
           any(start in thioester_matches for start in [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    ]
    
    if len(acyl_chain_bonds) != 1:
        return False, f"Found {len(acyl_chain_bonds)} double bonds in the acyl chain, need exactly 1"

    return True, "Structure is a monounsaturated fatty acyl-CoA with one double bond in the acyl chain"