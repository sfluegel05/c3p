"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:32847 dodecanoate ester
Defined as: Any fatty acid ester in which the carboxylic acid component is lauric acid (dodecanoic acid).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find dodecanoate fragment
    dodecanoate_fragment = Chem.MolFromSmiles("CCCCCCCCCCCC(=O)O")
    dodecanoate_matches = mol.GetSubstructMatches(dodecanoate_fragment)
    if not dodecanoate_matches:
        return False, "Dodecanoate fragment not found"

    # Check for ester bonds involving the dodecanoate fragment as the carboxylic acid component
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    dodecanoate_ester_matches = []
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[1])
        for dodecanoate_match in dodecanoate_matches:
            if ester_atom.GetIdx() in dodecanoate_match:
                dodecanoate_carbonyl_idx = [idx for idx in dodecanoate_match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8][0]
                dodecanoate_carbonyl_atom = mol.GetAtomWithIdx(dodecanoate_carbonyl_idx)
                if mol.GetBondBetweenAtoms(dodecanoate_carbonyl_idx, ester_atom.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    dodecanoate_ester_matches.append(match)
                    break

    if not dodecanoate_ester_matches:
        return False, "Dodecanoate fragment not involved in an ester bond as the carboxylic acid component"

    return True, "Contains dodecanoate fragment as the carboxylic acid component of an ester bond"