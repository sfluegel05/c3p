"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: CHEBI:35558 diketone
A diketone is a compound that contains two ketone functionalities.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdMolTransforms

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ketone groups
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Check for at least two ketone groups
    if len(ketone_matches) < 2:
        return False, f"Found {len(ketone_matches)} ketone groups, need at least 2"

    # Check if ketone groups are part of the same molecule
    # (not separate fragments)
    if not rdMolDescriptors.CalcNumAliphaticRings(mol) == 0:
        # Enumerates all possible tautomers
        tautomers = list(rdMolTransforms.EnumerateTautomers(mol))
        has_valid_tautomer = False
        for tautomer in tautomers:
            tautomer_matches = tautomer.GetSubstructMatches(ketone_pattern)
            if len(tautomer_matches) >= 2:
                has_valid_tautomer = True
                break
        if not has_valid_tautomer:
            return False, "Ketone groups are not part of the same molecule"

    # Additional checks to ensure the matched groups are indeed ketones
    # and not other functional groups like esters or carboxylic acids
    for atom_idx in set().union(*ketone_matches):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() != 0:
            return False, "Matched groups are not ketones"

    return True, "Contains at least two ketone functionalities"