"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    This class is characterized by the inclusion of a thioester-linked long-chain fatty acid with Coenzyme A.

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

    # Define Coenzyme A moiety pattern (simplified for illustrative purposes)
    coa_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Define thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Analyze carbon chain length to ensure it is between C13 and C22
    # (Assume C chain follows C=O in match)
    for match in thioester_matches:
        carbon_chain_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() > match[0]]
        carbon_count = sum(1 for atom_idx in carbon_chain_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
        if 13 <= carbon_count <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {carbon_count} C atoms in chain"

    return False, "Carbon chain length not in C13-C22 for long-chain fatty acid"