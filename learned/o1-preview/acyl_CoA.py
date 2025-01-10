"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester resulting from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the CoA moiety (excluding the thiol group)
    coa_smarts = Chem.MolFromSmarts("""
    NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O
    """)
    if coa_smarts is None:
        return False, "Error in defining CoA SMARTS pattern"

    # Define the SMARTS pattern for the thioester linkage (-C(=O)-S-)
    thioester_smarts = Chem.MolFromSmarts("C(=O)S")
    if thioester_smarts is None:
        return False, "Error in defining thioester SMARTS pattern"

    # Search for thioester linkage in the molecule
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "Thioester linkage not found"

    # Search for CoA moiety in the molecule
    if not mol.HasSubstructMatch(coa_smarts):
        return False, "Coenzyme A moiety not found"

    # Check if the thioester linkage connects the acyl group to CoA
    # Find the atoms involved in the thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    coa_matches = mol.GetSubstructMatches(coa_smarts)
    if not thioester_matches or not coa_matches:
        return False, "Thioester linkage or CoA moiety not properly matched"

    # Verify connectivity between thioester and CoA
    # Get the atom indices for sulfur in thioester and nitrogen in CoA
    sulfur_atoms = [atom_idx for match in thioester_matches for atom_idx in match if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'S']
    coa_atoms = [atom_idx for match in coa_matches for atom_idx in match]

    # Check if sulfur atom in thioester is connected to CoA moiety
    connected = False
    for s_idx in sulfur_atoms:
        for coa_idx in coa_atoms:
            path = Chem.rdmolops.GetShortestPath(mol, s_idx, coa_idx)
            if path and len(path) > 0:
                connected = True
                break
        if connected:
            break

    if not connected:
        return False, "Thioester linkage not connected to CoA moiety"

    return True, "Contains CoA moiety attached via thioester linkage to an acyl group"