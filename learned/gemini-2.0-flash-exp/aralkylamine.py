"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine where at least one of the alkyl substituents is an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for an amine (primary, secondary, or tertiary, but not quaternary)
    amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No amine group found"

    # 2. Check for at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("[c]")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
          return False, "No aromatic ring found"

    # 3. Check if any amine is connected to an aromatic ring via an alkyl chain.
    for amine_match in amine_matches:
        amine_atom_index = amine_match[0]
        amine_atom = mol.GetAtomWithIdx(amine_atom_index)

        # Function to recursively find connected atoms and check for aromatic carbon
        def find_connected_aromatic(current_atom, visited_atoms,depth):
            if current_atom.GetIdx() in visited_atoms:
                return False

            visited_atoms.add(current_atom.GetIdx())


            if current_atom.GetIsAromatic():
                 return True
            
            if depth > 10:  #Avoid infinite recursion on crazy molecules
                return False

            for neighbor in current_atom.GetNeighbors():
              if neighbor.GetAtomicNum() == 6:
                    if find_connected_aromatic(neighbor, visited_atoms,depth +1):
                        return True
            return False


        for neighbor in amine_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6: #Check carbons connected to the nitrogen
                if find_connected_aromatic(neighbor,set(),1):
                     return True, "Contains an amine group connected to an aromatic ring via an alkyl chain"
    
    return False, "No amine group connected to an aromatic ring via an alkyl chain"