"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
Definition: An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define the thioester linkage pattern (without CoA details)
    # Thioester linkage: C(=O)S
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not thioester_pattern:
        return False, "Invalid thioester SMARTS pattern"

    # Check for thioester linkage
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Define the 3-oxo-fatty acyl chain pattern
    # 3-oxo group on fatty acyl chain: [#6]-C(=O)-C-C(=O)-[#6]
    # This represents a carbon chain with a keto group at the 3-position
    three_oxo_acyl_smarts = "[#6]-C(=O)-C-C(=O)-[#6]"
    three_oxo_acyl_pattern = Chem.MolFromSmarts(three_oxo_acyl_smarts)
    if not three_oxo_acyl_pattern:
        return False, "Invalid 3-oxo-fatty acyl SMARTS pattern"

    # Check for 3-oxo-fatty acyl chain attached via thioester linkage
    matches = mol.GetSubstructMatches(three_oxo_acyl_pattern)
    if not matches:
        return False, "3-oxo-fatty acyl chain not found"

    # Check that the 3-oxo-fatty acyl chain is connected to the thioester linkage
    found_acyl_coa = False
    for match in matches:
        # Get the atom indices for the pattern match
        atom_indices = match
        # Check if one end is connected to the thioester carbonyl carbon
        acyl_carbon = atom_indices[1]  # The first carbonyl carbon in pattern
        # Get the atom corresponding to the acyl carbon
        acyl_atom = mol.GetAtomWithIdx(acyl_carbon)
        # Check if the acyl carbon is connected to a sulfur atom
        for neighbor in acyl_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'S':
                # Check if the sulfur is part of a thioester linkage
                bond = mol.GetBondBetweenAtoms(acyl_atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Found thioester linkage
                    found_acyl_coa = True
                    break
        if found_acyl_coa:
            break

    if not found_acyl_coa:
        return False, "3-oxo-fatty acyl chain not connected via thioester linkage"

    # Optionally, check the length of the fatty acyl chain
    # Count the number of carbons in the longest carbon chain
    fragments = Chem.GetMolFrags(mol, asMols=True)
    max_chain_length = 0
    for frag in fragments:
        atom_nums = [atom.GetAtomicNum() for atom in frag.GetAtoms()]
        c_count = atom_nums.count(6)
        if c_count > max_chain_length:
            max_chain_length = c_count
    if max_chain_length < 10:
        return False, f"Fatty acyl chain too short: {max_chain_length} carbons"

    return True, "Molecule is a 3-oxo-fatty acyl-CoA"