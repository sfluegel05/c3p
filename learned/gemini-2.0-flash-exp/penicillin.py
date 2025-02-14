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

    # 1. Penam core pattern (no specific numbering)
    penam_core_pattern = Chem.MolFromSmarts("[N]1[C]([C]2[S][C](C)(C)[C]2[C](O)=O)[C](=O)1")
    if penam_core_pattern is None:
        return False, "Invalid Penam Core SMARTS"
    
    substructure_matches = mol.GetSubstructMatches(penam_core_pattern)
    if not substructure_matches:
       return False, "Penam core not found"

    # 2. and 3. Two methyl substituents at position 2 and 3, and carboxylate at pos 3 are already captured by the SMARTS

    # 4. Carboxamido group at position 6 (R-CONH-).
    carboxamido_pattern = Chem.MolFromSmarts("[N][C](=O)")
    if carboxamido_pattern is None:
        return False, "Invalid carboxamido SMARTS"
    
    # Identify the carbon at position 6. Find the nitrogen of the core and then one of its neighbors should be the carbon we look for.
    match = substructure_matches[0]
    
    # The nitrogen in the core will be the first atom in the match.
    nitrogen_index = match[0]

    # Get the nitrogen of the core.
    nitrogen_atom = mol.GetAtomWithIdx(nitrogen_index)
    
    # The carbon at position 6 is a neighbor of the nitrogen in the core. There could be other neighbors, but it will always be a carbon.
    carbon_6_index = None
    for neighbor in nitrogen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            carbon_6_index = neighbor.GetIdx()
            break
    
    if carbon_6_index is None:
        return False, "Could not identify carbon at position 6"
    
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