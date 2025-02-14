"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
"""
from rdkit import Chem
from collections import deque

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine in which the alkyl group is substituted by an aromatic group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify amine nitrogens (exclude amides, nitriles, imines, etc.)
    amine_nitrogens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            is_amine = True
            # Exclude nitrogens with double or triple bonds (e.g., imines, amides, nitriles)
            for bond in atom.GetBonds():
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    is_amine = False
                    break
            if not is_amine:
                continue
            # Exclude quaternary ammonium (N with 4 bonds)
            if atom.GetTotalDegree() > 3:
                is_amine = False
            if is_amine:
                amine_nitrogens.append(atom)
    
    if not amine_nitrogens:
        return False, "No amine group found"
    
    # For each amine nitrogen, search for a path to an aromatic carbon via single bonds and carbons only
    for nitrogen in amine_nitrogens:
        visited = set()
        queue = deque()
        queue.append((nitrogen, 0))
        max_steps = 10  # Limit search to paths of max length 10
        found_aromatic = False
        while queue:
            current_atom, steps = queue.popleft()
            if steps > max_steps:
                continue
            if current_atom.GetIdx() in visited:
                continue
            visited.add(current_atom.GetIdx())
            # Skip the nitrogen atom itself
            if current_atom.GetIdx() != nitrogen.GetIdx():
                # Check if current atom is an aromatic carbon
                if current_atom.GetAtomicNum() == 6 and current_atom.GetIsAromatic():
                    found_aromatic = True
                    break
                # Only traverse carbons
                if current_atom.GetAtomicNum() != 6:
                    continue
            # Traverse neighbors via single bonds
            for bond in current_atom.GetBonds():
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                neighbor = bond.GetOtherAtom(current_atom)
                if neighbor.GetIdx() in visited:
                    continue
                queue.append((neighbor, steps+1))
        if found_aromatic:
            return True, "Amine group connected via alkyl chain to an aromatic ring"
    return False, "No aralkylamine substructure found"