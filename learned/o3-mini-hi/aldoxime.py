"""
Classifies: CHEBI:22307 aldoxime
"""
#!/usr/bin/env python
"""
Classifies: Oximes of aldehydes (aldoxime).
Definition: Aldoximes are the oximes derived from aldehydes, having the structure R–CH=N–OH,
where the carbon (CH) bears exactly one hydrogen (the aldehydic hydrogen).
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime (an oxime of an aldehyde) based on its SMILES string.
    
    In an aldoxime the functional group is R–CH=N–OH, meaning the carbon in the C=N bond
    must have exactly one hydrogen (from the original aldehyde) and the nitrogen is bound to an -OH group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains an aldoxime functional group, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are accurate.
    mol = Chem.AddHs(mol)
    
    # Define the SMARTS pattern for an aldoxime group.
    # The pattern [C;H1](=[N][O;H1]) looks for a carbon atom with exactly one hydrogen,
    # double bonded to a nitrogen atom; that nitrogen is bonded to an oxygen atom which itself has one hydrogen (OH group).
    aldoxime_pattern = Chem.MolFromSmarts("[C;H1](=[N][O;H1])")
    if aldoxime_pattern is None:
        return False, "Error creating SMARTS pattern for aldoxime"
    
    # Check for the aldoxime substructure in the molecule.
    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Molecule contains an aldoxime functional group (R-CH=N-OH)"
    else:
        return False, "Aldoxime functional group (R-CH=N-OH) not found"

# End of program.