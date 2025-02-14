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
    An 11,12-saturated fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) molecule in which the bond between
    carbons 11 and 12 in the fatty acyl chain is saturated (a single bond).

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

    # Use a more general CoA SMARTS pattern (excluding stereochemistry)
    coa_smarts = 'NC(=O)CCNC(=O)C(O)C(C)(C)CO[P](=O)(O)O[P](=O)(O)OCC1OC(CO[P](=O)(O)O)C(O)C1O'
    coa_mol = Chem.MolFromSmarts(coa_smarts)
    if coa_mol is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found"

    # Find the thioester linkage (C(=O)S)
    thioester_smarts = 'C(=O)S'
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_mol)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Identify the carbonyl carbon in the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]  # Index of the carbonyl carbon

    # Traverse the fatty acyl chain starting from the carbonyl carbon
    fatty_acyl_chain = [carbonyl_c_idx]
    current_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
    prev_atom_idx = None

    while True:
        neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom_idx]
        if not neighbors:
            break
        next_atom = neighbors[0]
        fatty_acyl_chain.append(next_atom.GetIdx())
        prev_atom_idx = current_atom.GetIdx()
        current_atom = next_atom

    # The fatty acyl chain should have at least 12 carbons to check the bond between carbons 11 and 12
    if len(fatty_acyl_chain) < 12:
        return False, "Fatty acyl chain is shorter than 12 carbons"

    # Get the bond between carbons 11 and 12
    c11_idx = fatty_acyl_chain[10]  # Zero-based indexing
    c12_idx = fatty_acyl_chain[11]
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