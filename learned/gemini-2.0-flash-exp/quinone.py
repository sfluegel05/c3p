"""
Classifies: CHEBI:36141 quinone
"""
ring_mol = Chem.Mol(mol)
        remove_atoms_ids = []
        for atom in mol.GetAtoms():
           if not atom.GetIdx() in ring:
             remove_atoms_ids.append(atom.GetIdx())

        ring_mol = Chem.RWMol(mol)
        remove_atoms_ids.sort(reverse=True)
        for atom_id in remove_atoms_ids:
            ring_mol.RemoveAtom(atom_id)

        ring_mol = ring_mol.GetMol()

        #check if the ring is aromatic
        if ring_mol.GetRingInfo().IsRingAromatic(ring):