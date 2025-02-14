"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nitrohydrocarbon, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for at least one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro group found"

    # 2. Verify that all nitrogens are part of a -NO2 group attached to a carbon
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 7: # Nitrogen
        is_nitro = False
        for match in nitro_matches:
          if atom.GetIdx() in match:
            is_nitro = True
            nitrogen_idx = atom.GetIdx()
            carbon_neighbor_count = 0;
            for neighbor in mol.GetAtomWithIdx(nitrogen_idx).GetNeighbors():
              if neighbor.GetAtomicNum() == 6:
                carbon_neighbor_count+=1
            if carbon_neighbor_count != 1:
              return False, "Nitro group is not bonded to exactly one carbon"
            break # already checked this nitrogen
        if not is_nitro:
            return False, "Nitrogen is not part of a -NO2 group"


    # 3. Check if the remaining structure is a hydrocarbon after removing nitro groups
    
    # This part is changed to a substructure approach. We first identify all atoms that are part of the -NO2 group.
    nitro_atoms = set()
    for match in nitro_matches:
      for idx in match:
         nitro_atoms.add(idx)
    
    # create a copy of the molecule to check if its a hydrocarbon after removing NO2 atoms
    mol_copy = Chem.Mol(mol)
    
    # We will make a copy and remove only the nitro atoms. Then check the result
    edit_mol = Chem.EditableMol(mol_copy)
    nitro_atoms = sorted(list(nitro_atoms), reverse=True) # important to remove in reversed order
    
    for atom_idx in nitro_atoms:
       edit_mol.RemoveAtom(atom_idx)
    mol_copy = edit_mol.GetMol()

    # Now check if all remaining atoms are either carbon or hydrogen
    for atom in mol_copy.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
          return False, "After removing nitro groups, non-hydrocarbon fragment remains"
    
    
    return True, "Molecule is a nitrohydrocarbon"