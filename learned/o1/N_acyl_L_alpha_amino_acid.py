"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is any L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol)
    Chem.rdCIPLabeler.AssignCIPLabels(mol)

    # Search for chiral carbons with 'S' configuration
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.HasProp('_CIPCode') and atom.GetProp('_CIPCode') == 'S':
            # Potential alpha carbon
            alpha_c = atom

            # Initialize flags
            has_carboxyl = False
            has_N_acyl = False

            # Get neighbors
            neighbors = alpha_c.GetNeighbors()
            for nbr in neighbors:
                if nbr.GetAtomicNum() == 7:
                    # Neighbor is nitrogen
                    nitrogen = nbr
                    # Check if nitrogen is acylated (connected to C=O)
                    is_acylated = False
                    for n_nbr in nitrogen.GetNeighbors():
                        if n_nbr.GetAtomicNum() == 6 and n_nbr != alpha_c:
                            # Check if carbon is a carbonyl carbon (C=O)
                            carbonyl_c = n_nbr
                            has_c_double_o = False
                            for bond in carbonyl_c.GetBonds():
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    atom1 = bond.GetBeginAtom()
                                    atom2 = bond.GetEndAtom()
                                    if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8) or
                                        (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6)):
                                        has_c_double_o = True
                            if has_c_double_o:
                                is_acylated = True
                                break
                    if is_acylated:
                        has_N_acyl = True
                elif nbr.GetAtomicNum() == 6:
                    # Check for carboxyl group attached to alpha carbon
                    carboxyl_c = nbr
                    # Check if this carbon is connected to two oxygens (carboxyl group)
                    o_count = 0
                    for cc_nbr in carboxyl_c.GetNeighbors():
                        if cc_nbr.GetAtomicNum() == 8:
                            o_count += 1
                    if o_count >= 2:
                        has_carboxyl = True

            if has_carboxyl and has_N_acyl:
                return True, "Molecule is an N-acyl-L-alpha-amino acid"

    return False, "Molecule is not an N-acyl-L-alpha-amino acid"