"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA results from the condensation of Coenzyme A (CoA) with a 3-hydroxy fatty acid via a thioester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simplified SMARTS pattern for the CoA moiety
    # Looking for the adenine ring connected to ribose and diphosphate
    coa_smarts = 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(O)(=O)OP(O)(=O)[O-])C(O)C3O'
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error parsing CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define a SMARTS pattern for the thioester linkage connecting CoA to acyl chain
    thioester_smarts = 'C(=O)SCCNC(=O)CCNC(=O)'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Error parsing thioester SMARTS pattern"

    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage connecting CoA to acyl chain found"

    # Define a SMARTS pattern for 3-hydroxy fatty acyl chain
    # Looking for a chain with a hydroxy group on the third carbon from the carbonyl carbon
    hydroxy_acyl_smarts = 'C(=O)C[C@H](O)C'
    hydroxy_acyl_pattern = Chem.MolFromSmarts(hydroxy_acyl_smarts)
    if hydroxy_acyl_pattern is None:
        return False, "Error parsing hydroxy acyl SMARTS pattern"

    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "No 3-hydroxy group on acyl chain found"

    # Optional: Check that the acyl chain is a fatty acid (long chain)
    # Count the number of carbons in the acyl chain
    # Here, we arbitrarily define a fatty acid as having at least 4 carbons in the chain
    carbonyl_c = mol.GetSubstructMatch(Chem.MolFromSmarts('C(=O)S'))[0]
    chain_lengths = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() == carbonyl_c or bond.GetEndAtomIdx() == carbonyl_c:
            neighbor_atom_idx = bond.GetOtherAtomIdx(carbonyl_c)
            chain_length = 0
            visited = set()
            stack = [neighbor_atom_idx]
            while stack:
                atom_idx = stack.pop()
                if atom_idx in visited:
                    continue
                visited.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    chain_length += 1
                    for nbr in atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx != carbonyl_c and nbr_idx not in visited:
                            stack.append(nbr_idx)
            chain_lengths.append(chain_length)
    if not any(length >= 4 for length in chain_lengths):
        return False, "Acyl chain is too short to be a fatty acid"

    return True, "Contains 3-hydroxy fatty acyl-CoA structure"