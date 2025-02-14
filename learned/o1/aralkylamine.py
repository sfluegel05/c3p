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
    
    # Identify amine nitrogens (include primary, secondary, tertiary, and quaternary amines)
    amine_nitrogens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            is_amine = True
            # Exclude nitrogens with double or triple bonds (e.g., amides, nitriles)
            for bond in atom.GetBonds():
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    is_amine = False
                    break
            if not is_amine:
                continue
            amine_nitrogens.append(atom)
    
    if not amine_nitrogens:
        return False, "No amine group found"
    
    # For each amine nitrogen, check for aralkylamine substructure
    for nitrogen in amine_nitrogens:
        # For each atom directly connected to nitrogen
        for bond in nitrogen.GetBonds():
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            neighbor = bond.GetOtherAtom(nitrogen)
            if neighbor.GetAtomicNum() != 6:
                continue  # Proceed only if neighbor is a carbon
            # Start path traversal from this carbon
            visited = set()
            queue = [(neighbor, 1)]  # Start with chain length 1
            max_chain_length = 3  # Limit the alkyl chain length
            found_aromatic = False
            while queue:
                current_atom, steps = queue.pop(0)
                if steps > max_chain_length:
                    continue
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                
                # Check if current atom is an aromatic carbon
                if current_atom.GetAtomicNum() == 6 and current_atom.GetIsAromatic():
                    found_aromatic = True
                    break
                # Only traverse carbons via single bonds
                for nbr_bond in current_atom.GetBonds():
                    if nbr_bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue
                    nbr_atom = nbr_bond.GetOtherAtom(current_atom)
                    if nbr_atom.GetAtomicNum() != 6:
                        continue
                    if nbr_atom.GetIdx() in visited:
                        continue
                    queue.append((nbr_atom, steps + 1))
            if found_aromatic:
                return True, "Amine group connected via alkyl chain to an aromatic ring"
    return False, "No aralkylamine substructure found"