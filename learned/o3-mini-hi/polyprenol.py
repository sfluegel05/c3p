"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: Polyprenol [Any member of the class of prenols possessing the general formula 
H-[CH2C(Me)=CHCH2]nOH, with n > 1 (more than one isoprene unit)].
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol must have a terminal primary alcohol (–CH2OH) and a carbon backbone composed 
    of at least two repeating isoprene units. The isoprene unit is represented as CH2–C(CH3)=CH–CH2.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple of classification result and reason.
    """
    
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule has at least one terminal primary alcohol.
    # Primary alcohol pattern: a CH2 group linked to an –OH
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "Molecule does not contain a terminal primary alcohol group"
    
    # Define a SMARTS pattern for a typical isoprene repeating unit:
    # (CH2)−C(CH3)=CH−CH2
    isoprene_unit = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]")
    
    # Look for substructure matches of the isoprene repeating unit.
    matches = mol.GetSubstructMatches(isoprene_unit)
    
    # Count the number of non-overlapping occurrences.
    # (A simple approach is to use the match count even if atoms are reused; 
    # polyprenols should contain a series of such fragments.)
    n_units = len(matches)
    
    # We require at least 2 isoprene units (the definition: composed of more than one isoprene unit).
    if n_units < 2:
        return False, f"Only {n_units} isoprene unit(s) detected; at least 2 are required"
    
    # Optionally check molecular weight to ensure size is compatible 
    # with typical polyprenols (expected to be in the upper range)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for a polyprenol"
    
    return True, f"Contains terminal primary alcohol and at least {n_units} isoprene units"

# For testing, one might call:
if __name__ == "__main__":
    # Examples of polyprenols (a few provided SMILES examples)
    test_smiles = [
        "CC(C)=CCC\\C(C)=C\\CO",  # geraniol (1 isoprene unit actually)
        "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",  # (2E,6E,10E)-omega-hydroxyfarnesol (3 isoprene units)
        "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO"  # (E,E,E)-geranylgeraniol (~4 isoprene units)
    ]
    for sm in test_smiles:
        res, reason = is_polyprenol(sm)
        print(f"SMILES: {sm}\n-> polyprenol: {res}, Reason: {reason}\n")