"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohol
Definition:
    A secondary alcohol is defined as a compound in which every –OH group 
    (after explicitly adding hydrogens) is attached to a saturated (sp³) carbon 
    that is bonded to exactly two other carbon atoms and one hydrogen.
    
    In addition, if disqualifying functional groups such as carboxylic acid/carboxylate 
    or phosphate groups (and certain carboxamide groups that are not directly linked to 
    the secondary alcohol carbon) are present, the molecule is not considered a pure 
    secondary alcohol.
    
Our approach:
  1. Parse the SMILES and add explicit hydrogens.
  2. Check for disqualifying substructures – for example a carboxylic acid (or its salt) 
     or a phosphate group. For carboxamides we allow only those directly adjacent to a secondary‐OH
     (e.g. in lactamide).
  3. Next, iterate over all oxygen atoms that may be part of an –OH group (i.e. atoms with exactly 
     two neighbors, one hydrogen and one carbon). For each, check that its attached carbon is sp³ 
     and (ignoring the –OH oxygen) has exactly two carbon neighbors and one hydrogen.
  4. If at least one –OH is encountered and every –OH conforms to the secondary pattern then 
     the molecule is classified as a secondary alcohol.
"""

from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    Here a “secondary alcohol” means that every hydroxyl group (–OH) in the molecule 
    is attached to a saturated (sp³) carbon that is bonded to exactly two other carbons
    and one hydrogen—and the molecule does not contain disqualifying functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a secondary alcohol, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can accurately count attached hydrogens.
    mol = Chem.AddHs(mol)
    
    # --- STEP 1: DISQUALIFY MOLECULES WITH UNWANTED FUNCTIONAL GROUPS ---
    # Carboxylic acids (or their deprotonated forms): e.g. ...C(=O)[OH] or ...C(=O)[O-]
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    if mol.HasSubstructMatch(acid_smarts):
        return False, "Contains carboxylic acid or carboxylate group"
    
    # Phosphate group: common pattern P(=O)(O)(O)
    phosphate_smarts = Chem.MolFromSmarts("P(=O)(O)(O)")
    if mol.HasSubstructMatch(phosphate_smarts):
        return False, "Contains phosphate group"
    
    # Carboxamide groups – allow only when the –OH is directly attached to the same carbon
    # (as in lactamide). Generic carboxamide pattern:
    carboxamide_smarts = Chem.MolFromSmarts("[CX3](=O)N")
    # Allowed pattern: a secondary alcohol carbon directly adjacent to a carboxamide;
    # e.g. CH(O)C(=O)N
    allowed_amide_smarts = Chem.MolFromSmarts("[C;X4]([#6])([#6])[H]O[C](=O)N")
    if mol.HasSubstructMatch(carboxamide_smarts) and (not mol.HasSubstructMatch(allowed_amide_smarts)):
        return False, "Contains carboxamide group not adjacent to a secondary alcohol"
    
    # --- STEP 2: IDENTIFY AND EVALUATE ALL -OH GROUPS ---
    total_OH = 0
    secondary_OH = 0
    
    # Iterate over all atoms looking for oxygen atoms that might be in an -OH group.
    for atom in mol.GetAtoms():
        # We only care about oxygen atoms.
        if atom.GetAtomicNum() != 8:
            continue
        
        # In an -OH group (after adding hydrogens), oxygen should have exactly 2 neighbors:
        # one hydrogen and one other atom.
        if atom.GetDegree() != 2:
            continue
        
        neighbors = atom.GetNeighbors()
        h_neighbor = None
        c_neighbor = None
        for nb in neighbors:
            if nb.GetAtomicNum() == 1:
                h_neighbor = nb
            elif nb.GetAtomicNum() == 6:
                c_neighbor = nb
        # If not exactly one hydrogen and one carbon are found, skip this oxygen.
        if h_neighbor is None or c_neighbor is None:
            continue
        
        # We have identified an -OH group.
        total_OH += 1
        
        # Check that the carbon attached to the –OH is sp³.
        if c_neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue  # This -OH is not on a saturated carbon.
        
        # Now, for a secondary alcohol the attached carbon (excluding the -OH oxygen)
        # should have exactly: 2 carbon neighbors and 1 hydrogen.
        count_carbon = 0
        count_hydrogen = 0
        for nb in c_neighbor.GetNeighbors():
            # Skip the oxygen that defines the -OH group.
            if nb.GetIdx() == atom.GetIdx():
                continue
            if nb.GetAtomicNum() == 6:
                count_carbon += 1
            elif nb.GetAtomicNum() == 1:
                count_hydrogen += 1
            # Other atoms (e.g. heteroatoms) make the substitution pattern deviate from the ideal.
        if count_carbon == 2 and count_hydrogen == 1:
            secondary_OH += 1
            
    if total_OH == 0:
        return False, "No identifiable -OH groups (in proper -OH form) found"
    
    if total_OH == secondary_OH:
        return True, f"All {total_OH} alcohol group(s) are secondary (each -OH is on a sp³ carbon bound to 2 carbons and 1 hydrogen)"
    else:
        return False, (f"Not all identified alcohol groups are secondary "
                       f"(found {secondary_OH} secondary out of {total_OH} total -OH groups)")

# Example usage:
# Uncomment the following lines to test a few examples:
# test_smiles = "CC(O)CC"   # Butan-2-ol, a classical secondary alcohol.
# print(is_secondary_alcohol(test_smiles))