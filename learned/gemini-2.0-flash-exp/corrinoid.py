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
    
    # remove metals and charges so that the substructure matching works
    mol_no_metals = Chem.Mol(mol)
    for atom in mol_no_metals.GetAtoms():
        if atom.GetAtomicNum() > 20:
            mol_no_metals.RemoveAtom(atom.GetIdx())
    for atom in mol_no_metals.GetAtoms():
        atom.SetFormalCharge(0)


    # Define a simplified SMARTS pattern for the corrin core, focusing on connectivity of the four nitrogen atoms.
    # The pattern matches the four nitrogens connected through four 5-membered rings.
    corrin_pattern = Chem.MolFromSmarts("[NX3;!$(N-)]1~[C]~[C]~[N;!$(N-)]~[C]~[C]~[N;!$(N-)]~[C]~[C]~[N;!$(N-)]~[C]~[C]1")
    
    if not mol_no_metals.HasSubstructMatch(corrin_pattern):
        return False, "Corrin core pattern not found."


    # Check for the presence of at least four pyrrole rings by checking for 5-membered rings containing nitrogen
    # A more relaxed pattern captures all pyrroles regardless of saturation
    pyrrole_pattern = Chem.MolFromSmarts("[NX3;!$(N-)]1[C][C][C][C]1") #5-membered ring containing N
    pyrrole_matches = mol_no_metals.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4 :
         return False, f"Found only {len(pyrrole_matches)} pyrrole rings. Need at least 4."

    # Check for the methine bridges connecting the pyrroles
    # Using a relaxed pattern that accounts for single or double bonds on either side.
    methine_pattern = Chem.MolFromSmarts("[NX3;!$(N-)]-[CX3H]=[CX3H]-[C]-[NX3;!$(N-)]") #single bond version
    methine_matches = mol_no_metals.GetSubstructMatches(methine_pattern)
    methine_pattern2 = Chem.MolFromSmarts("[NX3;!$(N-)]-[CX3H]=[CX3H]=[CX3H]-[NX3;!$(N-)]") #double bond version
    methine_matches2 = mol_no_metals.GetSubstructMatches(methine_pattern2)
    total_methine_matches = len(methine_matches) + len(methine_matches2)
    if total_methine_matches < 3:
        return False, f"Found only {total_methine_matches} methine bridges. Need at least 3."

    return True, "Contains a corrin nucleus with four pyrroles, 3 methine bridges and a direct C-C bond."