"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2, where R can be H.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific carboxamidine pattern [C](=[N])[N]
    core_pattern = Chem.MolFromSmarts("[C](=[N])[N]")
    matches = mol.GetSubstructMatches(core_pattern)
    
    if not matches:
        return False, "Core carboxamidine structure C(=N)-N not found."

    for match in matches:
        # Get the indices of the matched atoms
        c_idx, n1_idx, n2_idx = match

        # Check that carbon has exactly one double bond to N and 3 other single bonds 
        c_atom = mol.GetAtomWithIdx(c_idx)
        c_degree = c_atom.GetDegree()
        c_bonds = [bond for bond in c_atom.GetBonds() if bond.GetOtherAtomIdx(c_idx) != n1_idx] # ignore bond to N1
        c_valence = c_atom.GetTotalValence()
        if c_degree != 3:  # check explicit connections
            return False, f"Carbon {c_idx} does not have 3 single bonds."
        if c_valence > 4:
            return False, f"Carbon {c_idx} is bonded to more than 4 atoms"
        if c_bonds[0].GetBondType() == Chem.rdchem.BondType.DOUBLE:  # check no other double bonds
            return False, f"Carbon {c_idx} has an unexpected double bond"
        if c_bonds[1].GetBondType() == Chem.rdchem.BondType.DOUBLE:  # check no other double bonds
            return False, f"Carbon {c_idx} has an unexpected double bond"

        # Check that N1 is double bonded to C and singly bonded to something else
        n1_atom = mol.GetAtomWithIdx(n1_idx)
        n1_degree = n1_atom.GetDegree()
        n1_bonds = n1_atom.GetTotalValence()
        if n1_degree != 2:
            return False, f"Nitrogen {n1_idx} does not have 2 bonds"
        if n1_bonds != 3:
          return False, f"Nitrogen {n1_idx} is bonded to more than 3 atoms"
        
        
        # Check that N2 is connected to the C via the N
        n2_atom = mol.GetAtomWithIdx(n2_idx)
        n2_degree = n2_atom.GetDegree()
        n2_bonds = n2_atom.GetTotalValence()
        
        
        if n2_bonds > 3 : # only 2 atoms, could be 3 if N+, but we will not include N+ since it is part of the second pattern
           return False, f"Nitrogen {n2_idx} is bonded to more than 3 atoms"
        
        
        # Verify not Guanidine C(=N)(N)N
        guanidine_pattern = Chem.MolFromSmarts("[C](=[N])([N])[N]")
        if mol.HasSubstructMatch(guanidine_pattern):
           return False, "Molecule is a Guanidine"

    return True, "Molecule matches the core structure of a carboxamidine"