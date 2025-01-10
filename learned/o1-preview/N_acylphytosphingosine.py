"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:76955 N-acylphytosphingosine
"""

from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine
    having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acylphytosphingosine
    # Pattern for amide bond with nitrogen attached to phytosphingosine backbone
    n_acylphytosphingosine_smarts = Chem.MolFromSmarts("""
    [#6](=O)-[NX3;H0;!$(N-C=O)]-[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC
    """)

    if not n_acylphytosphingosine_smarts:
        return False, "Failed to create SMARTS pattern"

    # Search for matches
    matches = mol.GetSubstructMatches(n_acylphytosphingosine_smarts)
    if matches:
        # Verify long aliphatic chains
        for match in matches:
            # Indices of matched atoms
            carbonyl_c_idx = match[0]
            nitrogen_idx = match[1]
            c1_idx = match[2]  # C attached to nitrogen

            # Check acyl chain length (fatty acyl group)
            acyl_chain = Chem.rdmolops.GetShortestPath(mol, carbonyl_c_idx, nitrogen_idx)
            acyl_chain_length = len(acyl_chain) - 2  # Exclude C=O and N atoms
            if acyl_chain_length < 10:
                continue  # Short acyl chain, not typical for N-acylphytosphingosine

            # Check sphingoid base chain length
            c1_atom = mol.GetAtomWithIdx(c1_idx)
            sphingoid_chain_length = 0
            visited = set()
            atoms_to_visit = [n.GetIdx() for n in c1_atom.GetNeighbors() if n.GetIdx() != nitrogen_idx]
            while atoms_to_visit:
                current_idx = atoms_to_visit.pop()
                if current_idx in visited:
                    continue
                visited.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)
                if current_atom.GetAtomicNum() == 6:  # Carbon
                    sphingoid_chain_length += 1
                    neighbors = [n.GetIdx() for n in current_atom.GetNeighbors() if n.GetIdx() not in visited]
                    atoms_to_visit.extend(neighbors)
            if sphingoid_chain_length < 12:
                continue  # Short sphingoid base chain

            return True, "Molecule matches N-acylphytosphingosine pattern"

    return False, "Molecule does not match N-acylphytosphingosine pattern"