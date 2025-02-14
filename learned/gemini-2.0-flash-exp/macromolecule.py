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
    
    # Calculate the number of atoms in a chain (simplified version)
    def get_longest_chain_fast(mol):
       if not mol:
           return 0
       max_chain_len = 0

       for a1 in mol.GetAtoms():
          if a1.GetAtomicNum() != 1:
              for a2 in mol.GetAtoms():
                if a2.GetAtomicNum() != 1:
                  try:
                    path = Chem.GetShortestPath(mol,a1.GetIdx(),a2.GetIdx())
                    max_chain_len = max(max_chain_len,len(path))
                  except:
                     pass
       return max_chain_len

    num_atoms_in_longest_chain = get_longest_chain_fast(mol)
    
    
    # Apply heuristics for macromolecule classification
    if mol_wt > 1500:
        return True, "Very high molecular weight, likely a macromolecule."
    
    if num_heavy_atoms > 80 or num_rotatable_bonds > 25 or num_atoms_in_longest_chain > 30:
          return True, "High molecular weight or large number of heavy atoms or rotatable bonds or atoms in longest chain, indicating a macromolecule."
    
    elif num_heavy_atoms > 100 and num_rotatable_bonds > 30:
          return True, "Large number of heavy atoms and rotatable bonds, indicating a macromolecule"
    elif num_rotatable_bonds > 50 or num_atoms_in_longest_chain > 50:
           return True, "High number of rotatable bonds or atoms in chain, indicating a macromolecule."

    else:
        return False, "Does not meet macromolecule criteria based on molecular weight, heavy atoms, number of rotatable bonds or atoms in chain."