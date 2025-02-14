"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:64381 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine has a glycerol backbone with acyl groups at positions 1 and 2,
    and a phospho-L-serine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern (without stereochemistry to be more inclusive)
    glycerol_pattern = Chem.MolFromSmarts("C(C(O[*]))(C(O[*]))CO[P](=O)(O)O[C](C(N)C(=O)O)")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "Glycerol backbone with phosphate and serine not found"

    # Check for ester linkages attached to glycerol carbons
    ester_pattern = Chem.MolFromSmarts("C(OC(=O)[#6])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups attached to glycerol, need 2"

    # Check for phospho-L-serine group
    phospho_serine_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC[C](N)C(=O)O")
    if not mol.HasSubstructMatch(phospho_serine_pattern):
        return False, "Phospho-L-serine group not found"

    # Check that ester-linked chains are fatty acids (long hydrocarbon chains)
    acyl_chain_lengths = []
    ester_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O[C]"))
    for bond in ester_bonds:
        acyl_carbon_idx = bond[0]  # Carbonyl carbon
        # Trace the chain from the carbonyl carbon
        chain_length = 0
        visited = set()
        stack = [acyl_carbon_idx]
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                for neighbor in atom.GetNeighbors():
                    nbr_idx = neighbor.GetIdx()
                    nbr_atomic_num = neighbor.GetAtomicNum()
                    if nbr_idx not in visited and nbr_atomic_num in (6, 1):  # Continue through carbons and hydrogens
                        stack.append(nbr_idx)
        acyl_chain_lengths.append(chain_length)

    if len(acyl_chain_lengths) < 2:
        return False, "Less than two acyl chains found"
    if min(acyl_chain_lengths) < 8:
        return False, "Acyl chains are too short to be fatty acids"

    return True, "Molecule is a 3-sn-phosphatidyl-L-serine with correct glycerol backbone and substituents"