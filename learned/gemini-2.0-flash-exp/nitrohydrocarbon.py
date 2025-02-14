"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon with one or more hydrogen atoms replaced by a nitro group (-NO2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nitrohydrocarbon, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for invalid atoms. Only C, H, N, O are allowed
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 7, 8]: # 1(H), 6(C), 7(N), 8(O)
            return False, f"Contains invalid atom: {atom.GetSymbol()}"

    # 2. Check for at least one nitro group
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"

    # 3. Create a copy of the molecule, remove the nitro groups and check that what remains is a hydrocarbon
    mol_copy = Chem.Mol(mol) # create a copy of the molecule
    nitro_matches = mol_copy.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro group found"
    
    # remove the nitro groups from the copy
    for match in nitro_matches:
        # get all the atoms involved in the substructure match
        atoms_to_remove = [mol_copy.GetAtomWithIdx(idx).GetIdx() for idx in match] 
        
        # Get the atom that is not nitrogen or oxygen (the carbon)
        # Note: This is valid since we know that the nitrogen is always 
        # bonded to two oxygens and one carbon, and because the substructure 
        # match must be the -NO2 group.
        carbon_idx = -1;
        for atom_idx in match:
            atom = mol_copy.GetAtomWithIdx(atom_idx);
            if atom.GetSymbol() == 'C':
                carbon_idx = atom_idx;
                break;
        
        if carbon_idx == -1:
          return False, "Could not find carbon bonded to -NO2 group"
        
        bonds_to_remove = []
        # Find all the bonds we need to remove by comparing with original graph
        for bond in mol_copy.GetBonds():
           if bond.GetBeginAtomIdx() in atoms_to_remove and bond.GetEndAtomIdx() in atoms_to_remove:
               bonds_to_remove.append(bond.GetIdx())
           elif bond.GetBeginAtomIdx() == carbon_idx or bond.GetEndAtomIdx() == carbon_idx:
             begin_atom_idx = bond.GetBeginAtomIdx()
             end_atom_idx = bond.GetEndAtomIdx()
             if (begin_atom_idx in atoms_to_remove and end_atom_idx not in atoms_to_remove) or \
                (end_atom_idx in atoms_to_remove and begin_atom_idx not in atoms_to_remove):
                    bonds_to_remove.append(bond.GetIdx())
        
        bonds_to_remove.sort(reverse=True)
        for bond_idx in bonds_to_remove:
            mol_copy.RemoveBond(mol_copy.GetBondWithIdx(bond_idx).GetBeginAtomIdx(), mol_copy.GetBondWithIdx(bond_idx).GetEndAtomIdx())
        
        atoms_to_remove.sort(reverse=True)
        for atom_idx in atoms_to_remove:
            mol_copy.RemoveAtom(atom_idx)
    
    for atom in mol_copy.GetAtoms():
      if atom.GetAtomicNum() not in [1,6]:
        return False, "After removing nitro groups, non-hydrocarbon fragment remains"

    # 4. Check that all nitrogens are part of a -NO2 group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            if atom.GetTotalDegree() != 3: # must be bound to 3 other atoms
                return False, "Nitrogen is not part of -NO2 group (does not have 3 bonds)"
            
            oxygen_count = 0;
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    oxygen_count += 1
            if oxygen_count != 2: # must be bonded to two oxygens
                return False, "Nitrogen is not part of -NO2 group (does not have 2 oxygen neighbors)"

    return True, "Molecule is a nitrohydrocarbon"