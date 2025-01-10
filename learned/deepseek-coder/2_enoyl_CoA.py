"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: CHEBI:77557 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety (more flexible pattern)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for thioester bond in the acyl chain
    thioester_pattern = Chem.MolFromSmarts("[CX3]=[OX1]S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond found in the acyl chain"

    # Check for double bond between positions 2 and 3 in the acyl chain
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) == 0:
        return False, "No double bond between positions 2 and 3 in the acyl chain"

    # Ensure the double bond is between positions 2 and 3 relative to the thioester bond
    for thioester_match in thioester_matches:
        thioester_atom = thioester_match[1]  # The sulfur atom in the thioester bond
        for double_bond_match in double_bond_matches:
            atom1, atom2 = double_bond_match[0], double_bond_match[1]
            # Check if the double bond is between positions 2 and 3 relative to the thioester bond
            if mol.GetBondBetweenAtoms(atom1, thioester_atom) or mol.GetBondBetweenAtoms(atom2, thioester_atom):
                return True, "Contains CoA moiety with a double bond between positions 2 and 3 in the acyl chain attached via a thioester bond"

    return False, "No double bond between positions 2 and 3 in the acyl chain"