"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:25367 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a large molecule built from repeating smaller units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Check number of heavy (non-H) atoms
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    # Calculate number of rotatable bonds - indicates chain length
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Calculate the number of atoms in a chain
    
    def get_longest_chain(mol):
        """
         This function returns the number of heavy atoms in the longest chain in the molecule.
        """
        if not mol:
            return 0
        
        max_chain_length = 0
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() !=1 and a2.GetAtomicNum() != 1: # check that it is heavy atom bond.
              
              
              q = [(a1, [a1])]
              visited = set()
              while q:
                  curr_atom, path = q.pop(0)
                  
                  visited.add(curr_atom.GetIdx())
                  
                  max_chain_length = max(max_chain_length, len(path))
                  
                  for neighbor in curr_atom.GetNeighbors():
                     if neighbor.GetAtomicNum() != 1 and neighbor.GetIdx() not in visited:
                        q.append((neighbor, path + [neighbor]))
        return max_chain_length    
    num_atoms_in_longest_chain = get_longest_chain(mol)
    
    
    # Apply heuristics for macromolecule classification
    if mol_wt > 1000 or num_heavy_atoms > 80 or num_rotatable_bonds > 25 or num_atoms_in_longest_chain > 30:
          return True, "High molecular weight or large number of heavy atoms or rotatable bonds or atoms in longest chain, indicating a macromolecule."
    elif mol_wt > 1500:
          return True, "Very high molecular weight, likely a macromolecule."
    elif num_heavy_atoms > 100 and num_rotatable_bonds > 30:
          return True, "Large number of heavy atoms and rotatable bonds, indicating a macromolecule"
    elif num_rotatable_bonds > 50 or num_atoms_in_longest_chain > 50:
           return True, "High number of rotatable bonds or atoms in chain, indicating a macromolecule."

    else:
        return False, "Does not meet macromolecule criteria based on molecular weight, heavy atoms, number of rotatable bonds or atoms in chain."