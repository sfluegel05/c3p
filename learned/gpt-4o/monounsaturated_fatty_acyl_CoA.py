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
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for Coenzyme A (CoA) pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A pattern found"

    # Thioester pattern to find the start of the acyl chain
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found, implying absence of acyl chain"

    acyl_start_atom_idx = thioester_matches[0][0]  # The C involved in C=O

    # Use RDKit functionality to ascertain the acyl chain
    acyl_chain_end = False
    num_double_bonds = 0

    atom = mol.GetAtomWithIdx(acyl_start_atom_idx)
    visited_atoms = set([acyl_start_atom_idx])

    def traverse_chain(atom, visited):
        nonlocal num_double_bonds, acyl_chain_end
        
        # Traverse the atom's neighbors
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            idx = neighbor.GetIdx()
            if idx in visited:
                continue

            # Detect double bonds
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                num_double_bonds += 1

            visited.add(idx)

            # A simple heuristic to end traversal when CoA structure might get revisited
            if neighbor.GetAtomicNum() in (8, 7):  # Oxygen or nitrogen
                acyl_chain_end = True
                break
            
            traverse_chain(neighbor, visited)
            if acyl_chain_end:
                break

    traverse_chain(atom, visited_atoms)

    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds in the acyl chain, need exactly 1"

    return True, "Structure is a monounsaturated fatty acyl-CoA with one double bond in the acyl chain"