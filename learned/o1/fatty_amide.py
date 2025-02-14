"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:24038 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define amide functional group pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    
    # Find all amide groups in the molecule
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amides = len(amide_matches)
    
    if num_amides == 0:
        return False, "No amide functional group found"
    
    # For each amide group, check for an acyl chain with sufficient length
    for amide_match in amide_matches:
        carbonyl_c_idx = amide_match[0]
        o_idx = amide_match[1]
        n_idx = amide_match[2]
        
        # Exclude carbonyl oxygen and amide nitrogen from traversal
        exclude_atoms = set([o_idx, n_idx])
        
        # Get the acyl chain attached to the carbonyl carbon via carbon-carbon bonds
        chain_atoms = set()
        atoms_to_visit = [carbonyl_c_idx]
        visited_atoms = set()
        
        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                chain_atoms.add(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if nbr_idx in exclude_atoms or nbr_idx in visited_atoms:
                    continue
                if nbr.GetAtomicNum() == 6:
                    atoms_to_visit.append(nbr_idx)
        
        chain_length = len(chain_atoms)
        if chain_length >= 8:
            return True, f"Molecule is a fatty amide with acyl chain length {chain_length}"
    
    return False, "No acyl chain of sufficient length found"

# Example usage:
# smiles = "CCCCCCCCCCCCCCCC(=O)NCCO"  # N-palmitoyl ethanolamide
# result, reason = is_fatty_amide(smiles)
# print(result, reason)