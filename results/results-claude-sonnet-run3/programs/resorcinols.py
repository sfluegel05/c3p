from rdkit import Chem
from rdkit.Chem import AllChem

def is_resorcinols(smiles: str):
    """
    Determines if a molecule is a resorcinol (benzenediol with hydroxy groups meta to each other).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a resorcinol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for resorcinol (1,3-dihydroxybenzene)
    # Matches a benzene ring with two OH groups in meta position
    resorcinol_pattern = Chem.MolFromSmarts('c1c(O)cc(O)cc1')
    
    if mol.HasSubstructMatch(resorcinol_pattern):
        # Find all matches
        matches = mol.GetSubstructMatches(resorcinol_pattern)
        
        # Verify that the matched atoms form a proper resorcinol structure
        for match in matches:
            # Get the matched atoms
            ring_atoms = set(match)
            oh_positions = []
            
            # Find the OH groups
            for atom_idx in ring_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                        oh_positions.append(atom_idx)
            
            # If we found exactly two OH groups
            if len(oh_positions) == 2:
                # Calculate distance between OH groups in the ring
                path_len = len(Chem.GetShortestPath(mol, oh_positions[0], oh_positions[1]))
                # Meta position = 4 bonds apart in ring
                if path_len == 4:
                    substituents = []
                    for atom_idx in ring_atoms:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIdx() not in ring_atoms and neighbor.GetSymbol() != 'O':
                                substituents.append(neighbor.GetSymbol())
                    
                    if substituents:
                        return True, f"Substituted resorcinol with substituents: {', '.join(set(substituents))}"
                    else:
                        return True, "Unsubstituted resorcinol"
    
    return False, "No resorcinol pattern found"
# Pr=None
# Recall=None