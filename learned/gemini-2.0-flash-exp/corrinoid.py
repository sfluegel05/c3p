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

    # Remove charges before removing metals
    mol_no_charges = Chem.Mol(mol)
    for atom in mol_no_charges.GetAtoms():
        atom.SetFormalCharge(0)


    # Create an EditableMol to remove metal atoms
    mol_editable = Chem.EditableMol(mol_no_charges)
    atoms_to_remove = []

    for atom in mol_editable.GetMol().GetAtoms():
            if atom.GetAtomicNum() > 20:
                atoms_to_remove.append(atom.GetIdx())
    
    for idx in reversed(atoms_to_remove): #remove atoms in reverse so the indices do not get messed up when removing several atoms
        mol_editable.RemoveAtom(idx)
    
    mol_no_metals = mol_editable.GetMol()



    # Define a simplified SMARTS pattern for the corrin core, focusing on connectivity of the four nitrogen atoms.
    corrin_pattern = Chem.MolFromSmarts("[N]~[C]~[C]~[N]~[C]~[C]~[N]~[C]~[C]~[N]")
    if not mol_no_metals.HasSubstructMatch(corrin_pattern):
        return False, "Corrin core pattern not found."


    # Check for the presence of at least four pyrrole rings by checking for 5-membered rings containing nitrogen
    pyrrole_pattern = Chem.MolFromSmarts("[N]1[C][C][C][C]1") #5-membered ring containing N
    pyrrole_matches = mol_no_metals.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) < 4 :
         return False, f"Found only {len(pyrrole_matches)} pyrrole rings. Need at least 4."

    # Check for the methine bridges connecting the pyrroles. Using a relaxed pattern
    methine_pattern = Chem.MolFromSmarts("[N]-[C]=[C]-[C]-[N]")
    methine_matches = mol_no_metals.GetSubstructMatches(methine_pattern)
    if len(methine_matches) < 3:
         return False, f"Found only {len(methine_matches)} methine bridges. Need at least 3."


    return True, "Contains a corrin nucleus with four pyrroles, 3 methine bridges and a direct C-C bond."