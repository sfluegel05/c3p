"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine has a nitrogen atom bonded to two hydrogens and one sp3-hybridized carbon (NH2 group attached to a hydrocarbyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to accurately count hydrogen atoms
    mol = Chem.AddHs(mol)
    
    found_primary_amine = False
    
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen
        if atom.GetAtomicNum() == 7:
            # Check if nitrogen has degree 1 (connected to one heavy atom)
            if atom.GetDegree() == 1:
                # Check if nitrogen has two hydrogens
                if atom.GetTotalNumHs(includeNeighbors=True) == 2:
                    neighbor = atom.GetNeighbors()[0]
                    # Check if the neighbor atom is carbon
                    if neighbor.GetAtomicNum() == 6:
                        # Check if the bond between N and C is single
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Check if carbon is sp3-hybridized (tetrahedral)
                            if neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                                # Exclude if nitrogen is in a ring
                                if not atom.IsInRing():
                                    # Exclude if carbon is connected to any heteroatoms other than N or H (e.g., carbonyl group)
                                    carbon_heteroatoms = [nbr.GetAtomicNum() for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() not in (1, 6, 7)]
                                    if not carbon_heteroatoms:
                                        # Exclude if carbon is double-bonded to any atom (e.g., C=O)
                                        has_double_bond = False
                                        for nbr in neighbor.GetNeighbors():
                                            bond_to_nbr = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx())
                                            if bond_to_nbr.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                                has_double_bond = True
                                                break
                                        if not has_double_bond:
                                            # Exclude if carbon is aromatic
                                            if not neighbor.GetIsAromatic():
                                                found_primary_amine = True
                                                return True, "Contains primary amine group (NH2 attached to hydrocarbyl group)"
    if found_primary_amine:
        return True, "Contains primary amine group (NH2 attached to hydrocarbyl group)"
    else:
        return False, "No primary amine group found"

__metadata__ = {
    'chemical_class': {
        'name': 'primary amine',
        'definition': 'A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.'
    },
    'config': {
        'additional_info': 'Classifies primary amines based on NH2 group attached to a non-aromatic sp3 hybridized carbon, excluding amides and other functional groups.'
    }
}