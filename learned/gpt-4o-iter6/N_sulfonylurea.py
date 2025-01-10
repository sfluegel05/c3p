"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea contains a urea group where one of the nitrogen hydrogens
    is replaced by a sulfonyl group (-S(=O)(=O)-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
  
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    # Urea group: N-C(=O)-N (simplified depiction for pattern matching)
    urea_pattern = Chem.MolFromSmarts("NC(=O)N")  
    # N-sulfonyl group: N-S(=O)(=O)
    sulfonyl_pattern = Chem.MolFromSmarts("[NX3,NX4][SX4](=[OX1])(=[OX1])")

    # Check urea group presence
    if not mol.HasSubstructMatch(urea_pattern):
        return False, "No urea group found"
    
    # Check sulfonyl group substitution on nitrogen of urea
    if mol.HasSubstructMatch(sulfonyl_pattern):
        # Ensure the sulfonyl group is attached to one of the nitrogens of the urea
        for match in mol.GetSubstructMatches(urea_pattern):
            # Check the nitrogens of the matched urea structure
            urea_nitrogen_atoms = [match[0], match[2]]
            for n_idx in urea_nitrogen_atoms:
                ref_atom = mol.GetAtomWithIdx(n_idx)
                for neighbor in ref_atom.GetNeighbors():
                    if all([neighbor.GetSymbol() == "S", 
                            any(bond.GetBondTypeAsDouble() == 2.0 for bond in neighbor.GetBonds())]):
                        return True, "Contains N-sulfonylurea moiety"

    return False, "Does not contain an N-sulfonylurea moiety"