"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is characterized by a central alpha carbon connected to a protonated amino group and a deprotonated carboxylate group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found = False
    for atom in mol.GetAtoms():
        # Look for positively charged nitrogen atoms [NH3+]
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Get the carbon atoms connected to this nitrogen
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    alpha_carbon = neighbor
                    # Check if alpha carbon is connected to a carboxylate group
                    has_carboxylate = False
                    for neighbor2 in alpha_carbon.GetNeighbors():
                        if neighbor2.GetAtomicNum() == 6 and neighbor2.GetIdx() != atom.GetIdx():
                            # Potential carboxylate carbon
                            carboxylate_carbon = neighbor2
                            oxygens = [n for n in carboxylate_carbon.GetNeighbors()
                                       if n.GetAtomicNum() == 8]
                            if len(oxygens) == 2:
                                charges = [o.GetFormalCharge() for o in oxygens]
                                bonds = [carboxylate_carbon.GetBondBetweenAtoms(carboxylate_carbon.GetIdx(), o.GetIdx()).GetBondType()
                                         for o in oxygens]
                                if (set(charges) == {-1, 0} and
                                    Chem.rdchem.BondType.DOUBLE in bonds and
                                    Chem.rdchem.BondType.SINGLE in bonds):
                                    has_carboxylate = True
                                    break
                    if has_carboxylate:
                        found = True
                        break
            if found:
                break

    if found:
        return True, "Contains alpha-amino-acid zwitterion motif"
    else:
        return False, "Missing alpha-amino-acid zwitterion motif"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57421',
                          'name': 'alpha-amino-acid zwitterion',
                          'definition': 'An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group of any alpha-amino acid; major species at pH 7.3.',
                          'parents': ['CHEBI:58245', 'CHEBI:59871']},
    'message': None,
    'success': None}