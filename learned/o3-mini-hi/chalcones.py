"""
Classifies: CHEBI:23086 chalcones
"""
#!/usr/bin/env python3
"""
Classifies: Chalcones – compounds containing a 1,3-diphenylpropenone (benzylideneacetophenone) core or its dihydro derivatives.
That is, molecules containing an unsaturated ketone motif (Ar–C(=O)–CH=CH–Ar) 
or the dihydro variant (Ar–C(=O)–CH2–CH2–Ar), with both aromatic termini being benzene rings.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule belongs to the chalcone class (or chalcone derivative)
    based on its SMILES string. Our approach is to search for one of two core substructures:
       (i) an α,β-unsaturated ketone – benzylideneacetophenone -- with connectivity:
           benzene ring - C(=O) - CH=CH - benzene ring
       (ii) a dihydrochalcone variant with connectivity:
           benzene ring - C(=O) - CH2 - CH2 - benzene ring
    Both terminal aromatic groups are required to be a benzene ring.
    This method is less restrictive than checking if the carbonyl and linker atoms are non-cyclic.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a chalcone (or chalcone derivative), False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the chalcone cores.
    # (i) Unsaturated chalcone: aromatic ring - carbonyl - alkene - aromatic ring.
    unsat_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)C=CC2=ccccc2")
    # (ii) Dihydrochalcone: aromatic ring - carbonyl - CH2 - CH2 - aromatic ring.
    sat_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)CCc2ccccc2")
    
    # Check if the molecule matches either core.
    if mol.HasSubstructMatch(unsat_pattern):
        return True, "Contains chalcone core (α,β-unsaturated ketone variant with terminal benzene rings)"
    if mol.HasSubstructMatch(sat_pattern):
        return True, "Contains chalcone core (dihydrochalcone variant with terminal benzene rings)"
    
    # If no match was found, return a negative result.
    return False, "Chalcone core not found (expected connectivity: benzene ring - C(=O) - vinyl or ethylene linker - benzene ring)"

# For basic testing you could uncomment the lines below:
# test_smiles = "O=C(\\C=C\\c1ccccc1)c1ccccc1"  # trans-chalcone
# result, reason = is_chalcones(test_smiles)
# print(result, reason)