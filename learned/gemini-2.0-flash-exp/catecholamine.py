"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a benzene ring with hydroxyl groups at positions 1 and 2 (a catechol),
    and a 2-aminoethyl group attached at position 4 or 5, with possible substitutions.

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

    # Check for catechol moiety (1,2-dihydroxybenzene)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(*)c(*)c1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety (1,2-dihydroxybenzene) found"
    
    # Check for 2-aminoethyl group attached to the benzene ring, excluding positions 1 and 2
    #The substructure that fits is: benzene ring with 1,2-OH, and -CH2-CH2-NH2 attached to it
    full_pattern = Chem.MolFromSmarts("c1c(O)c(O)c([CH2][CH2][NH2])c(*)c1;!@c1c(O)c([CH2][CH2][NH2])c(*)c(O)c1;!@c1c(O)c(*)c([CH2][CH2][NH2])c(O)c1")
    
    if not mol.HasSubstructMatch(full_pattern):
        return False, "No 2-aminoethyl group correctly attached to the catechol ring found"

    return True, "Contains a catechol moiety and a 2-aminoethyl group attached to the ring"