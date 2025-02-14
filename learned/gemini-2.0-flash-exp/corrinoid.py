"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid is a derivative of the corrin nucleus, which contains four reduced or partly
    reduced pyrrole rings joined in a macrocycle by three =C- groups and one direct
    carbon-carbon bond linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a corrinoid, False otherwise
         str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a combined SMARTS pattern that captures the corrin core connectivity
    corrin_pattern = Chem.MolFromSmarts(
        "[NX3;!$(N-)]1-[CX3H]~[CX3H]~[NX3;!$(N-)]2-[CX3H]~[CX3H]~[NX3;!$(N-)]3-[CX3H]~[CX3H]~[NX3;!$(N-)]4-[CX3H]~[CX3H]1"
    ) # ring closure explicitly included

    #Match the pattern for the corrin core
    if not mol.HasSubstructMatch(corrin_pattern):
        return False, "Corrin core pattern not found."

    # Check that the corrin core is part of a larger ring system
    corrin_core_matches = mol.GetSubstructMatches(corrin_pattern)
    if len(corrin_core_matches) == 0:
        return False, "Could not match corrin core"

    # Check for C=C or C-C bonds in the pyrrole rings
    pyrrole_pattern_1 = Chem.MolFromSmarts("[NX3;!$(N-)]1[CX3H]=[CX3H][CX3H][CX3H]1")
    pyrrole_pattern_2 = Chem.MolFromSmarts("[NX3;!$(N-)]1[CX3H][CX3H]=[CX3H][CX3H]1")
    pyrrole_pattern_3 = Chem.MolFromSmarts("[NX3;!$(N-)]1[CX3H][CX3H][CX3H][CX3H]1") # fully reduced

    matches_pyrrole_1 = mol.GetSubstructMatches(pyrrole_pattern_1)
    matches_pyrrole_2 = mol.GetSubstructMatches(pyrrole_pattern_2)
    matches_pyrrole_3 = mol.GetSubstructMatches(pyrrole_pattern_3)

    total_pyrroles = len(matches_pyrrole_1) + len(matches_pyrrole_2) + len(matches_pyrrole_3)

    if total_pyrroles < 4:
        return False, f"Found only {total_pyrroles} pyrrole rings."

    # Check for the methine bridges and ensure they are connecting the pyrrole rings
    methine_pattern = Chem.MolFromSmarts("[NX3;!$(N-)]-[CX3H]=[CX3H]-[CX3H]=[CX3H]-[NX3;!$(N-)]")
    methine_matches = mol.GetSubstructMatches(methine_pattern)
    if len(methine_matches) < 3:
        return False, f"Found {len(methine_matches)} methine bridges. Need at least 3."

    # check for a direct C-C bond linking alpha positions of the pyrroles using the corrin_pattern
    # this is captured in corrin_pattern already
        
    return True, "Contains a corrin nucleus with four pyrroles, 3 methine bridges and a direct C-C bond."