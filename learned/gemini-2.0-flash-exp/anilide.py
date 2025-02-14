"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is an aromatic amide obtained by acylation of aniline.

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

    #  Find an amide bond
    amide_pattern = Chem.MolFromSmarts("[N]-[CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Find an aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("[c]")
    aromatic_matches = mol.GetSubstructMatches(aromatic_pattern)
    if not aromatic_matches:
        return False, "No aromatic ring found"


    # Check if the aromatic ring is directly or indirectly connected to the amide nitrogen
    connected = False
    for amide_match in amide_matches:
         for aromatic_match in aromatic_matches:
            for atom_idx in amide_match:
                amide_atom = mol.GetAtomWithIdx(atom_idx)
                if amide_atom.GetSymbol() == 'N':
                    for neighbor in amide_atom.GetNeighbors():
                         if neighbor.GetIdx() in aromatic_match:
                              connected = True
                              break
                         for n_neighbor in neighbor.GetNeighbors():
                             if n_neighbor.GetIdx() in aromatic_match:
                                  connected = True
                                  break
                    if connected:
                        break
            if connected:
                break


    if not connected:
         return False, "No direct or indirect connection between aromatic ring and amide nitrogen"

    return True, "Correct Anilide core detected"