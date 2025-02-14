"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid has a 3-ketone group and a C=C double bond between position 4 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for steroid core with explicit ring fusion and numbering (to check for positions)
    # Note that this SMARTS pattern specifically captures the fused ring system of steroids and does not
    # allow for any other ring system. It uses a more constrained definition of the core structure.
    steroid_core_pattern = Chem.MolFromSmarts("[C]12[C]([C])([C])[C]([C])([C])[C]3[C]1[C]([C])([C])[C]([C])([C])[C]4[C]2[C]([C])([C])[C]([C])([C])[C]34")
    
    # Define SMARTS for the 3-oxo and delta-4 double bond with explicit atoms positions
    #The 3-oxo-delta-4 pattern looks like C1-C(=O)-C=C-C, and the atoms positions correspond to atoms in the core_pattern
    oxo_group_pattern = Chem.MolFromSmarts("[C]1-[C](=[O])-[C]=[C]-[C]1")
    
    #Find if the substructures are present in the molecule.
    if not mol.HasSubstructMatch(steroid_core_pattern):
         return False, "No steroid core found"

    # We verify if the oxo group is present
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo-Delta(4) group found"

    # Find all the matches of the oxo group
    oxo_matches = mol.GetSubstructMatches(oxo_group_pattern)
    
    #Verify the location of the double bond and the carbonyl group.
    is_correct_location = False
    for match in oxo_matches:
        c_1 = match[0] # C that is part of the six membered ring, has to be at position 2
        c_2 = match[1] # Carbonyl carbon, has to be at position 3
        c_3 = match[2] # Double bond carbon, has to be at position 4
        c_4 = match[3] # Double bond carbon, has to be at position 5
        c_5 = match[4] # Carbon adjacent to double bond, has to be at position 6
        #Get the ring number of the atoms from the steroid core.
        for core_match in mol.GetSubstructMatches(steroid_core_pattern):
            core_c_2 = core_match[1]
            core_c_3 = core_match[2]
            core_c_4 = core_match[3]
            core_c_5 = core_match[4]
            
            if c_1 == core_c_2 and c_2 == core_c_3 and c_3 == core_c_4 and c_4 == core_c_5:
                is_correct_location = True
                break
        if is_correct_location:
          break #if the location is correct for one match, then the steroid is classified.

    if not is_correct_location:
        return False, "3-oxo-Delta(4) group not in correct location in steroid core."

    return True, "Contains steroid core with 3-oxo and Delta(4) double bond."