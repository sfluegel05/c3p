"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC)
Definition: “Any organic compound having an initial boiling point less than or equal 
to 250 °C (482 °F) measured at a standard atmospheric pressure of 101.3 kPa.”
Note: RDKit does not provide a direct method to predict boiling points. Accurately 
predicting boiling point generally requires a QSPR model or experimental data.
Thus, this implementation currently only verifies that the SMILES string is valid 
and represents an organic compound. The function then returns (None, None) 
to indicate that no rigorous classification has been performed.
"""
from rdkit import Chem

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its 
    estimated boiling point. According to the definition, a VOC must have an 
    initial boiling point <= 250 °C. However, as RDKit does not include boiling point 
    prediction, a proper classification cannot be performed.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool or None: True if the compound is classified as a VOC, False if not, 
                      or None if the prediction could not be performed.
        str or None: A reason for the classification (or lack thereof).
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic check: does the molecule contain carbon? (i.e. is it organic)
    if not any(atom.GetSymbol() == "C" for atom in mol.GetAtoms()):
        return False, "Not an organic compound (contains no carbon)"

    # Here we would predict the initial boiling point. Unfortunately, RDKit does not
    # include a boiling point estimator and reliable prediction requires a QSPR model.
    # Without such a model, we cannot decide whether the compound's boiling point is 
    # <= 250°C.
    
    return None, "Boiling point prediction not implemented; cannot classify VOC without QSPR model."