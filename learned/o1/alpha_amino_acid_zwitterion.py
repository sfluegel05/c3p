"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:57844 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is an amino acid zwitterion obtained by transfer of a proton from 
    the carboxy group to the amino group of any alpha-amino acid; major species at physiological pH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and add hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)

    alpha_amino_acid_found = False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # Skip if not carbon
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue  # Skip if not sp3 hybridized carbon

        neighbors = atom.GetNeighbors()
        nitrogen_neighbor = None
        carboxyl_neighbor = None

        for nbr in neighbors:
            if nbr.GetAtomicNum() == 7:
                # Found nitrogen neighbor
                nitrogen_neighbor = nbr
            elif nbr.GetAtomicNum() == 6:
                # Potential carboxyl carbon neighbor
                oxygens = [o for o in nbr.GetNeighbors() if o.GetAtomicNum() == 8]
                if len(oxygens) >= 2:
                    # Check for one double-bonded oxygen and one single-bonded oxygen
                    bond_types = [nbr.GetBondBetweenAtoms(nbr.GetIdx(), o.GetIdx()).GetBondType() for o in oxygens]
                    if Chem.rdchem.BondType.DOUBLE in bond_types and Chem.rdchem.BondType.SINGLE in bond_types:
                        carboxyl_neighbor = nbr

        if nitrogen_neighbor and carboxyl_neighbor:
            # Check that nitrogen has at least one hydrogen (implicit or explicit)
            n_hydrogens = nitrogen_neighbor.GetTotalNumHs()
            if n_hydrogens >= 1:
                alpha_amino_acid_found = True
                break

    if alpha_amino_acid_found:
        return True, "Contains alpha-amino-acid zwitterion structure"
    else:
        return False, "Alpha-amino-acid zwitterion pattern not found"