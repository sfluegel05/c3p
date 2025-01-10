"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation 
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Find phosphorus atoms in molecule
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    
    if not p_atoms:
        return False, "No phosphorus atoms found, not a phosphate"
        
    # For each phosphorus atom, check if it's part of a phosphate ester group
    for p_atom in p_atoms:
        # Check if P atom has 4 oxygen neighbors
        neighbors = p_atom.GetNeighbors()
        o_neighbors = [atom for atom in neighbors if atom.GetAtomicNum() == 8]
        if len(o_neighbors) != 4:
            continue  # Not a phosphate group
        # Check for one double bond (P=O)
        num_double_bonds = 0
        ester_o_atoms = []
        for o_atom in o_neighbors:
            bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), o_atom.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                num_double_bonds += 1
            else:
                # Single-bonded oxygen, possible ester oxygen
                # Check if connected to carbon
                o_neigh = [a for a in o_atom.GetNeighbors() if a.GetIdx() != p_atom.GetIdx()]
                for atom in o_neigh:
                    if atom.GetAtomicNum() == 6:
                        ester_o_atoms.append(o_atom)
        if num_double_bonds !=1:
            continue  # Not a typical phosphate group
        if not ester_o_atoms:
            continue  # No ester oxygen connected to carbon
        # For each ester oxygen atom connected to carbon, traverse the chain
        for ester_o_atom in ester_o_atoms:
            # Get the carbon atom connected to the ester oxygen
            o_neighbors = ester_o_atom.GetNeighbors()
            for neighbor in o_neighbors:
                if neighbor.GetIdx() == p_atom.GetIdx():
                    continue  # Skip P atom
                if neighbor.GetAtomicNum() == 6:
                    starting_c_atom = neighbor
                    break
            # Now traverse the chain starting from starting_c_atom
            visited_atoms = set()
            atoms_to_visit = [(starting_c_atom, None)]  # (current_atom, previous_atom)
            num_carbons = 0
            num_double_bonds_chain = 0
            num_methyl_branches = 0
            chain_atoms = []
            while atoms_to_visit:
                current_atom, previous_atom = atoms_to_visit.pop()
                if current_atom.GetIdx() in visited_atoms:
                    continue
                visited_atoms.add(current_atom.GetIdx())
                if current_atom.GetAtomicNum() != 6:
                    continue  # Skip non-carbon atoms
                num_carbons +=1
                chain_atoms.append(current_atom.GetIdx())
                # Check for double bonds between current atom and previous atom
                if previous_atom is not None:
                    bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), previous_atom.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        num_double_bonds_chain +=1
                # Check for methyl branches (carbon with degree 1)
                neighbors = current_atom.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != previous_atom.GetIdx():
                        if neighbor.GetDegree() == 1:
                            num_methyl_branches +=1
                        else:
                            # Continue traversal in linear fashion
                            atoms_to_visit.append((neighbor, current_atom))
            # Now, check the criteria for polyprenol chain
            if num_carbons >=20 and num_double_bonds_chain >=4 and num_methyl_branches >=4:
                if num_carbons %5 ==0:
                    return True, f"Contains phosphate ester group connected to a polyprenol chain of length {num_carbons} carbons"
                else:
                    return True, f"Contains phosphate ester group connected to a long chain of {num_carbons} carbons (not multiple of 5)"
            else:
                continue  # Try next phosphate group or ester oxygen
    return False, "No suitable phosphate ester group connected to polyprenol chain found"