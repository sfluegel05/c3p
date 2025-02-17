"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: Quaternary ammonium ion
Definition: A derivative of ammonium (NH4+) in which all four of the hydrogens bonded to nitrogen 
have been replaced with univalent (usually organyl) groups.
"""

from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a quaternary ammonium group based on its SMILES string.
    
    The function looks for at least one nitrogen atom carrying a positive charge with exactly
    four substituents (i.e. no attached hydrogens), which matches the definition of a quaternary ammonium ion.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains a quaternary ammonium group, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Option 1: Use a SMARTS pattern for a quaternary ammonium group.
    # The SMARTS "[N+;X4]" will match any nitrogen with a positive charge and 4 connections.
    pattern = Chem.MolFromSmarts("[N+;X4]")
    
    if mol.HasSubstructMatch(pattern):
        return True, "Found a nitrogen atom with +1 charge and four substituents (quaternary ammonium ion)"
    
    # Option 2: Alternatively we can iterate over atoms and check manually.
    for atom in mol.GetAtoms():
        # Check that the atom is a nitrogen with a +1 formal charge.
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Get the number of explicit neighbors. In a quaternary ammonium ion, there should be 4 bonds.
            if atom.GetDegree() == 4:
                # Ensure that no hydrogen is attached (quaternary means all hydrogens replaced).
                # Although RDKit typically includes implicit hydrogens if not explicitly drawn,
                # a formally charged quaternary ammonium group should not have any.
                if atom.GetTotalNumHs() == 0:
                    return True, "Found a nitrogen atom with +1 charge, four substituents, and no hydrogens (quaternary ammonium ion)"
    
    return False, "No quaternary ammonium group (N+ with 4 substituents) found"