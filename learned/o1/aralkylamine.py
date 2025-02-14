"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
"""
from rdkit import Chem

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
    
    # Identify amine nitrogens (exclude amides, nitriles, etc.)
    amine_nitrogens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            is_amine = True
            # Exclude nitrogens in amides (N-C(=O))
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    bonds = [bond for bond in nbr.GetBonds() if bond.GetOtherAtom(nbr).GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
                    if bonds:
                        is_amine = False
                        break
            # Exclude nitriles (C#N)
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                    is_amine = False
                    break
            if is_amine:
                amine_nitrogens.append(atom)
    
    if not amine_nitrogens:
        return False, "No amine group found"

    # For each amine nitrogen, search for a path to an aromatic ring via single bonds
    for nitrogen in amine_nitrogens:
        visited = set()
        stack = [(nitrogen, 0)]  # (current_atom, depth)
        found_aromatic = False
        max_depth = 50  # Prevent infinite loops, adjust if necessary
        
        while stack:
            current_atom, depth = stack.pop()
            if depth > max_depth:
                continue
            if current_atom.GetIdx() in visited:
                continue
            visited.add(current_atom.GetIdx())
            
            # Skip starting nitrogen for aromatic check
            if current_atom.GetIdx() != nitrogen.GetIdx():
                # Check if current atom is an aromatic carbon
                if current_atom.GetAtomicNum() == 6 and current_atom.GetIsAromatic():
                    found_aromatic = True
                    break
                # Exclude heteroatoms other than oxygen and sulfur in the chain
                elif current_atom.GetAtomicNum() not in [6, 8, 16]:
                    continue

            # Traverse neighbors via single bonds
            for bond in current_atom.GetBonds():
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue
                neighbor = bond.GetOtherAtom(current_atom)
                if neighbor.GetIdx() in visited:
                    continue
                # Avoid going back to the nitrogen if already passed
                if neighbor.GetAtomicNum() == 7 and neighbor != nitrogen:
                    continue
                stack.append((neighbor, depth + 1))
        if found_aromatic:
            return True, "Amine group connected via alkyl chain to an aromatic ring"

    return False, "No aralkylamine substructure found"