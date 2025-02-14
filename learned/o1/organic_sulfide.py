"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:25698 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide is defined as R-S-R, where R is not H (sulfur atom bonded to two carbon atoms via single bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over sulfur atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Atomic number of sulfur is 16
            # Check if sulfur has exactly two neighbors (degree 2)
            neighbor_atoms = atom.GetNeighbors()
            if len(neighbor_atoms) == 2:
                # Check if both neighbors are carbons
                if all(neighbor.GetAtomicNum() == 6 for neighbor in neighbor_atoms):
                    # Check that the bonds are single bonds
                    bonds = [mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()) for neighbor in neighbor_atoms]
                    if all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in bonds):
                        # Check that sulfur has no formal charge and a total valence of 2
                        if atom.GetFormalCharge() == 0 and atom.GetTotalValence() == 2:
                            return True, "Contains sulfur atom bonded to two carbons via single bonds (organic sulfide)"

    return False, "No sulfur atom bonded to two carbons via single bonds found"