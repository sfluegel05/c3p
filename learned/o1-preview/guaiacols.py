"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Iterate over atoms to find aromatic carbons with OH group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            # Check if this carbon is bonded to OH
            has_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        # Check if O atom has only one neighbor (it's an OH group)
                        if nbr.GetDegree() == 1:
                            has_OH = True
                            break
            if has_OH:
                # Now check for adjacent aromatic carbon with OCH3 group
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                        # Exclude the oxygen
                        if neighbor.GetIdx() != nbr.GetIdx():
                            for neigh_nbr in neighbor.GetNeighbors():
                                if neigh_nbr.GetAtomicNum() == 8:
                                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), neigh_nbr.GetIdx())
                                    if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                        # Check if O atom is connected to a methyl group
                                        if neigh_nbr.GetDegree() == 2:
                                            attached_carbons = [a for a in neigh_nbr.GetNeighbors() if a.GetIdx() != neighbor.GetIdx()]
                                            if len(attached_carbons) == 1:
                                                c_atom = attached_carbons[0]
                                                if c_atom.GetAtomicNum() == 6 and c_atom.GetDegree() == 1:
                                                    # Found OCH3 group at ortho position
                                                    return True, "Contains phenol with ortho-methoxy group"
    return False, "Does not contain phenol with ortho-methoxy group"