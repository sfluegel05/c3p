"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
"""
Classifies: 11,12-saturated fatty acyl-CoA(4-)
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    This means the molecule is a fatty acyl-CoA(4-) where the bond between carbons 11 and 12 
    in the fatty acyl chain is saturated (a single bond).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to accurately count atoms
    mol = Chem.AddHs(mol)

    # Use a comprehensive SMARTS pattern for the CoA moiety
    coa_smarts = '[NX3][CX3](=O)[CX4][CX4][NX3][CX3](=O)[CX4]([OX2H])[CX4]([CX4])([CX4])[CX4][OX2H][PX4](=O)([O-])[OX2][PX4](=O)([O-])[OX2][CX4][CX4]1[OX2][CX4]([OX2H])[CX4]1[OX2][PX4](=O)([O-])[O-]'
    coa_mol = Chem.MolFromSmarts(coa_smarts)
    if coa_mol is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found"

    # Find the thioester sulfur atom connected to the fatty acyl chain
    thioester_s_smarts = '[SX2][CX3](=O)[CX4]'
    thioester_s_mol = Chem.MolFromSmarts(thioester_s_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_s_mol)
    if not thioester_matches:
        return False, "Thioester linkage to CoA not found"

    # Identify the starting carbon of the fatty acyl chain
    fatty_acyl_start_idx = [match[2] for match in thioester_matches][0]

    # Traverse the fatty acyl chain starting from the alpha carbon
    fatty_acyl_atoms = []
    current_idx = fatty_acyl_start_idx
    visited = set()

    while True:
        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() != 6:
            break
        fatty_acyl_atoms.append(current_idx)
        visited.add(current_idx)
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                     if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]
        if not neighbors:
            break
        current_idx = neighbors[0]

    # The fatty acyl chain should have at least 12 carbons
    if len(fatty_acyl_atoms) < 12:
        return False, "Fatty acyl chain is shorter than 12 carbons"

    # Get the bond between carbons 11 and 12
    c11_idx = fatty_acyl_atoms[10]  # Zero-based indexing
    c12_idx = fatty_acyl_atoms[11]
    bond = mol.GetBondBetweenAtoms(c11_idx, c12_idx)
    if bond is None:
        return False, "No bond between carbons 11 and 12"

    bond_type = bond.GetBondType()
    if bond_type == Chem.rdchem.BondType.SINGLE:
        return True, "Bond between carbons 11 and 12 is saturated (single bond)"
    elif bond_type == Chem.rdchem.BondType.DOUBLE:
        return False, "Bond between carbons 11 and 12 is unsaturated (double bond)"
    else:
        return False, f"Bond between carbons 11 and 12 is neither single nor double bond (bond type: {bond_type})"