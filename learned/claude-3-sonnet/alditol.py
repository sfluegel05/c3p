"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:18237 alditol
An alditol is a carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH
(formally derivable from an aldose by reduction of the carbonyl group).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for acyclic structure
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1CCCCC1")):
        return False, "Molecule contains a ring system"

    # Check for terminal groups
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):
        return False, "Missing terminal -CH2OH group"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("OC")):
        return False, "Missing terminal -HOCH2 group"

    # Count carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check for -CH(OH)- pattern
    choh_pattern = Chem.MolFromSmarts("[CH](O)")
    choh_matches = mol.GetSubstructMatches(choh_pattern)
    n_choh = len(choh_matches)

    if c_count != n_choh + 2 or o_count != n_choh + 2:
        return False, "Incorrect number of carbon and oxygen atoms for alditol"

    # Check for linear arrangement of -CH(OH)- groups
    for idx in choh_matches:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetDegree() != 2:
            return False, "-CH(OH)- groups not linearly arranged"

    return True, "Acyclic polyol with correct terminal groups and linear -CH(OH)- arrangement"