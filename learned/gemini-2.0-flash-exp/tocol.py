"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton substituted at position 2 by a
    saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the chroman-6-ol core.
    # This pattern looks for a benzene ring fused to a dihydrofuran ring (ether), and a hydroxyl group at position 6.
    # The 6-membered ring is connected to the 2 position through a carbon atom.
    chromanol_pattern = Chem.MolFromSmarts('[c]1[c]([OX2])[c]([c])[c]2[O][CH2][CH2][c]2[c]1')
    core_matches = mol.GetSubstructMatches(chromanol_pattern)
    if not core_matches:
        return False, "Chroman-6-ol core not found"
    
    # Get the 2-position of the core.
    chromanol_match = core_matches[0]
    
    #Get the position 2 carbon
    position_2_carbon = None
    for atom_idx in chromanol_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == "C" and atom.GetExplicitValence() == 4:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "C" and neighbor.GetExplicitValence() == 4:
                     for next_neighbor in neighbor.GetNeighbors():
                        if next_neighbor.GetSymbol() == 'O':
                            position_2_carbon = atom
                            break
        if position_2_carbon is not None:
            break

    if position_2_carbon is None:
          return False, "Could not identify carbon at position 2."
    
    # Check for the isoprenoid chain and its length
    isoprenoid_chain_carbons = []
    visited_atoms = set()

    def traverse_chain(current_atom, chain):
        if current_atom.GetIdx() in visited_atoms:
            return
        
        visited_atoms.add(current_atom.GetIdx())
        if current_atom.GetSymbol() == "C":
          chain.append(current_atom)

        
        neighbors = current_atom.GetNeighbors()

        for neighbor in neighbors:
           if neighbor.GetSymbol() == "C" and neighbor.GetIdx() != position_2_carbon.GetIdx():
               traverse_chain(neighbor, chain)


    traverse_chain(position_2_carbon, isoprenoid_chain_carbons)

    isoprenoid_count = len(isoprenoid_chain_carbons)

    if isoprenoid_count != 15:
        return False, f"Isoprenoid chain does not have 15 carbons. Found {isoprenoid_count} carbon atoms in the chain."


    #Check the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300: # Tocols are generally above 300 Da
        return False, "Molecular weight too low for tocol"

    return True, "Contains chroman-6-ol core with a 15-carbon isoprenoid chain at position 2."