"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at the 16th carbon of a steroid core,
    where this hydroxy group has a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid core SMARTS pattern (more flexible and handles implicit Hs)
    steroid_core_smarts = "[C]1[C][C]2[C]3[C]([C]1)[C][C]4[C]3[C]([C]2)[C][C]5[C]4[C]([C])([C]5)"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    # Find substructure match for steroid core
    core_match = mol.GetSubstructMatch(steroid_core_pattern)

    if not core_match:
         return False, "No steroid core found"

    # Map the match to the steroid core
    core_atoms = {}
    for i, atom_idx in enumerate(core_match):
       core_atoms[i+1] = atom_idx  # 1-based indexing for core atoms
    
    # Identify C16 in the core
    # In the steroid pattern:
    # carbon 1: index 0
    # carbon 2: index 1
    # carbon 3: index 2
    # carbon 4: index 9
    # carbon 5: index 8
    # carbon 6: index 7
    # carbon 7: index 6
    # carbon 8: index 5
    # carbon 9: index 3
    # carbon 10: index 10
    # carbon 11: index 11
    # carbon 12: index 12
    # carbon 13: index 4
    # carbon 14: index 13
    # carbon 15: index 14
    # carbon 16: index 15
    # C16 is the 16th atom of the core, we get its index
    c16_idx = core_atoms.get(16)
    if c16_idx is None:
        return False, "Cannot determine index of C16"

    # Get the C16 atom object
    c16_atom = mol.GetAtomWithIdx(c16_idx)
    
    # Check for a hydroxyl group directly attached to C16 with beta configuration
    has_beta_oh = False

    for neighbor in c16_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8: # Oxygen atom
           o_atom = neighbor
           
           # Count H neighbors on the O atom. We expect 1 hydrogen neighbor
           h_count = 0
           for h_neighbor in o_atom.GetNeighbors():
               if h_neighbor.GetAtomicNum() == 1: # Hydrogen
                   h_count += 1
           if h_count == 1:
               o_c16_bond = mol.GetBondBetweenAtoms(c16_atom.GetIdx(), o_atom.GetIdx())
               if o_c16_bond:
                    if o_c16_bond.GetStereo() == Chem.BondStereo.STEREOUP: # check for beta config
                           has_beta_oh = True
                           break # found a suitable OH group, exit loop


    if not has_beta_oh:
        return False, "No 16-beta hydroxyl group found"
    
    return True, "16beta-hydroxy steroid found"