"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA (C13 to C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("[N]1C=NC2=C1N=CN=C2N[C@H]3[C@@H]([C@H]([C@@H](O3)COP(=O)(O)OP(=O)(O)OCC(O)C([NH])=O)O)OP(=O)(O)O")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No CoA substructure found"
    
    # Look for thioester bond
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"
    
    # Check for long-chain fatty acid (C13-C22)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]~[CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No long-chain fatty acid found"
    
    # Count double bonds in fatty acid chain
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtomIdx() in fatty_acid_matches[0] and bond.GetEndAtomIdx() in fatty_acid_matches[0])
    
    # Long-chain fatty acyl-CoAs typically have double bonds in the fatty acid chain
    if double_bonds == 0:
        return False, "No double bonds found in fatty acid chain"
    
    return True, "Contains CoA substructure with a long-chain fatty acid (C13-C22) attached via a thioester bond"