"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is a 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a general benzopyran core pattern
    # This pattern checks for a fused 6-membered aromatic ring (c1ccccc1) and a 6-membered ring containing an oxygen. 
    # X is used to denote any atom for now, C is aromatic carbon and the second bond of the oxygen atom should also
    # be to an aromatic carbon. Also, consider cases where there is no carbonyl.
    core_pattern1 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(-[*])-C=C2') # For 1-benzopyran
    core_pattern2 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(=O)-C-C=C2') # For 1-benzopyran-2-one (coumarin) 
    core_pattern3 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(-[*])-CC2') # For dihydro-1-benzopyran
    
    has_core = mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2) or mol.HasSubstructMatch(core_pattern3)
    if not has_core:
        return False, "Does not contain a 1-benzopyran core."


    # Now, check for an aryl group directly connected to the 4 position of the benzopyran.
    # The * matches any atom which is attached to the carbon of the core.
    aryl_substituent_pattern1 = Chem.MolFromSmarts('c1ccccc1-[*]') # A simple phenyl
    aryl_substituent_pattern2 = Chem.MolFromSmarts('c1ccccc1c1-[*]') # Biphenyl
    
    
    
    # Find the carbon at position 4 in the core.
    matches = []
    if mol.HasSubstructMatch(core_pattern1):
       matches = mol.GetSubstructMatches(core_pattern1)
       match_pattern = core_pattern1
    elif mol.HasSubstructMatch(core_pattern2):
       matches = mol.GetSubstructMatches(core_pattern2)
       match_pattern = core_pattern2
    elif mol.HasSubstructMatch(core_pattern3):
        matches = mol.GetSubstructMatches(core_pattern3)
        match_pattern = core_pattern3
    
    if not matches:
        return False, "Could not find core match"
   
    found_aryl_at_pos4 = False
    for match in matches:
        # the position of the carbon attached to the core is always at position 4
        # get the atom that corresponds to the [*] in the SMARTS. This will be the 
        # atom attached to C at position 4 of the benzopyran core.
        pos4_atom_idx = match_pattern.GetSubstructMatch(Chem.MolFromSmarts('*'))[0] #Atom index of substituent at pos4
        pos4_atom = mol.GetAtomWithIdx(match[pos4_atom_idx])
        
        # Now check neighbors and if any of them are attached to an aromatic ring
        for neighbor in pos4_atom.GetNeighbors():
           sub_mol = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, atomBonds = list(mol.GetBonds()), rootedAtAtom = neighbor.GetIdx(), singleBondsOnly = True))
           if sub_mol.HasSubstructMatch(aryl_substituent_pattern1) or sub_mol.HasSubstructMatch(aryl_substituent_pattern2):
                found_aryl_at_pos4 = True
                break
        if found_aryl_at_pos4:
            break #found the aryl group at position 4

    if not found_aryl_at_pos4:
         return False, "Does not have an aryl substituent at position 4 of the 1-benzopyran core."
        
    return True, "Contains a 1-benzopyran core with an aryl substituent at position 4."