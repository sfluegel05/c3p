"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid or a corrinoid precursor/derivative based on its SMILES string.
    Corrinoids are derivatives of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings
    joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid or corrinoid precursor/derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for corrin macrocycle pattern
    corrin_pattern = Chem.MolFromSmarts("[N&R]1[C&R]=2[N&R]=3[C&R]=4[N&R]=5[C&R]=6[N&R]=1[C&R]=2[C&R]=3[C&R]=4[C&R]=5[C&R]=6")
    corrin_matches = mol.GetSubstructMatches(corrin_pattern)

    # Look for precorrin and cobalt-precorrin patterns
    precorrin_pattern = Chem.MolFromSmarts("[N&R]1[C&R]=2[N&R]=3[C&R]=4[N&R]=5[C&R]=6[N&R]=1[C&R]=2[C&R]=3[C&R]=4[C&R]=5[C&R]=6[C&R]")
    precorrin_matches = mol.GetSubstructMatches(precorrin_pattern)

    # Look for cobalt coordination
    cobalt_pattern = Chem.MolFromSmarts("[Co]")
    cobalt_matches = mol.GetSubstructMatches(cobalt_pattern)

    # Check for corrin macrocycle or precursor/derivative pattern
    if corrin_matches or precorrin_matches:
        if cobalt_matches:
            return True, "Contains corrin macrocycle or precursor/derivative with cobalt coordination"
        else:
            return True, "Contains corrin macrocycle or precursor/derivative without cobalt coordination"

    # Check for common corrinoid modifications
    modifications = ["[N&R]C(=O)CC[C@@H]", "[N&R]C[C@@H](OP([O-])([O-])=O)", "[N&R]C[C@@H](OP(O)(O)=O)O[C@@H]", "[N&R]c1nc[nH]c1=O"]
    for mod in modifications:
        mod_pattern = Chem.MolFromSmarts(mod)
        mod_matches = mol.GetSubstructMatches(mod_pattern)
        if mod_matches:
            return True, f"Contains corrinoid modification pattern: {mod}"

    return False, "No corrinoid or corrinoid precursor/derivative pattern found"