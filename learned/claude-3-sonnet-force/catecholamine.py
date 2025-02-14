"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18148 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a 4-(2-aminoethyl)pyrocatechol (4-(2-aminoethyl)benzene-1,2-diol)
    or a derivative formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pyrocatechol core with 2-aminoethyl side chain
    catecholamine_pattern = Chem.MolFromSmarts("c1c(O)cc(O)c(c1)CCNC")
    catecholamine_matches = mol.GetSubstructMatches(catecholamine_pattern)
    if not catecholamine_matches:
        return False, "No catecholamine core found"

    # Check for substitutions on the pyrocatechol core
    substituted = False
    for match in catecholamine_matches:
        core_atom = mol.GetAtomWithIdx(match[0])
        if core_atom.GetDegree() > 3:
            substituted = True
            break

    # Check for additional hydroxyl groups
    num_hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if num_hydroxyl_groups > 2:
        substituted = True

    # Check for other substitutions on the catecholamine core
    catecholamine_core = Chem.MolFromSmarts("c1c(O)cc(O)c(CCNC)c1")
    other_substitutions = mol.HasSubstructMatch(catecholamine_core, useChirality=True)
    if other_substitutions:
        substituted = True

    if substituted:
        return True, "Substituted catecholamine derivative"
    else:
        return True, "Contains unsubstituted catecholamine core"