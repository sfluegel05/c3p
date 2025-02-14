"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:51443 anilide

An anilide is any aromatic amide obtained by acylation of aniline.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aniline substructure (aromatic ring with NH2 group)
    aniline_pattern = Chem.MolFromSmarts("c1ccc(cc1)N")
    aniline_matches = mol.GetSubstructMatches(aniline_pattern)

    # Look for aromatic amide group (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("C(=O)Nc")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    # Check if aniline and amide groups are connected
    connected = False
    for aniline_match in aniline_matches:
        for amide_match in amide_matches:
            aniline_atom = aniline_match[-1]  # Last atom in the aniline match
            amide_atom = amide_match[2]  # Nitrogen atom in the amide match
            if mol.GetBondBetweenAtoms(aniline_atom, amide_atom) is not None:
                connected = True
                break

    # Check if the aniline ring is aromatic
    aniline_ring = mol.GetAtomWithIdx(aniline_match[0]).GetIsAromatic()

    # Check if the amide group is part of an aromatic system
    amide_ring = False
    for atom in amide_match:
        if mol.GetAtomWithIdx(atom).GetIsAromatic():
            amide_ring = True
            break

    # Check for acylation pattern (C(=O)N-Ar)
    acylation_pattern = Chem.MolFromSmarts("C(=O)Nc")
    acylation_matches = mol.GetSubstructMatches(acylation_pattern)

    if connected and aniline_ring and amide_ring and acylation_matches:
        return True, "Contains an aromatic aniline group connected to an aromatic amide group via acylation"
    else:
        return False, "Does not match the structural requirements for an anilide"