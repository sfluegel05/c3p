"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one carbon-chlorine (C-Cl) bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carbon-chlorine bond (any bond type)
    c_cl_bond = Chem.MolFromSmarts("[#6][Cl]")
    if mol.HasSubstructMatch(c_cl_bond):
        # Ensure chlorine is not attached to non-carbon atoms
        matches = mol.GetSubstructMatches(c_cl_bond)
        for match in matches:
            atom1 = mol.GetAtomWithIdx(match[0])  # Carbon atom
            atom2 = mol.GetAtomWithIdx(match[1])  # Chlorine atom
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 17:
                return True, "Contains at least one carbon-chlorine (C-Cl) bond"
        return False, "No carbon-chlorine (C-Cl) bonds found"
    else:
        return False, "No carbon-chlorine (C-Cl) bonds found"