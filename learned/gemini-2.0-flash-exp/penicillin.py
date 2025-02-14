"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as having a penam core, two methyl substituents at position 2,
    a carboxylate substituent at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Penam core pattern (with numbering for easy referencing)
    #     7    6     5
    #   --N--[C@H]--[CH]--
    #  |   |     |     |
    #  |   C------C------S
    #  |   2     3     4
    #   --    C=O  --
    #
    penam_core_pattern = Chem.MolFromSmarts("[N1][C@H]([C@H]2[S][C@](C)(C)[C@H]2C(O)=O)C(=O)1")
    if penam_core_pattern is None:
        return False, "Invalid Penam Core SMARTS"
    
    substructure_matches = mol.GetSubstructMatches(penam_core_pattern)
    if not substructure_matches:
       return False, "Penam core not found"


    # 2. Two methyl substituents at position 2 and 3. Carboxylate substituent at position 3. These are already checked in the penam core definition above.
    # 4. Carboxamido group at position 6 (R-CONH-). It's not a specific R group.
    carboxamido_pattern = Chem.MolFromSmarts("[N][C](=O)")
    if carboxamido_pattern is None:
        return False, "Invalid carboxamido SMARTS"

    
    # get the penam substructure match, and from the match get the index of the atom at the 6 position of the penam core.
    match = substructure_matches[0]
    
    # get the index of the nitrogen of the core at position 7 (index 0 in the pattern)
    nitrogen_7_index = match[0]
    
    # get the index of the carbon at position 6 (index 1 in the pattern)
    carbon_6_index = match[1]

    # Check if any carboxamido group matches the whole molecule.
    all_carboxamido_matches = mol.GetSubstructMatches(carboxamido_pattern)

    if not all_carboxamido_matches:
        return False, "Carboxamido group not found"
    
    found_correct_carboxamido = False
    for carboxamido_match in all_carboxamido_matches:
        #check that any of the carboxamido nitrogen atoms are neighbours of the carbon 6 of the penam core.
        for atom_index in carboxamido_match:
           atom = mol.GetAtomWithIdx(atom_index)
           for neighbor in atom.GetNeighbors():
               if neighbor.GetIdx() == carbon_6_index:
                   found_correct_carboxamido = True
                   break
           if found_correct_carboxamido:
               break
        if found_correct_carboxamido:
            break
            
    if not found_correct_carboxamido:
        return False, "Carboxamido group not at position 6."
    

    return True, "Molecule contains penam core, two methyl substituents, carboxylate, and carboxamido group"