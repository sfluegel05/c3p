"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:24024 Monoterpenoid 
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is a terpenoid derived from a monoterpene, which has a C10 skeleton 
    that may be rearranged or modified by the removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count heavy atoms (non-hydrogen)
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    
    # Check if the molecule contains 10 heavy atoms (C10 skeleton)
    if n_heavy_atoms != 10:
        return False, f"Number of heavy atoms ({n_heavy_atoms}) does not match C10 skeleton"
    
    # Check for presence of carbon ring
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No carbon ring found"
    
    # Count carbon atoms
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 10:
        return False, f"Too few carbon atoms ({n_carbons}) for monoterpenoid"
    
    # Check for unsaturated bonds
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    n_triple_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.TRIPLE)
    n_unsaturated_bonds = n_double_bonds + n_triple_bonds
    if n_unsaturated_bonds > 2:
        return False, "Too many unsaturated bonds for monoterpenoid"
    
    # Check for functional groups (alcohols, esters, ethers, aldehydes, ketones, carboxylic acids)
    for atom in mol.GetAtoms():
        atom_num = atom.GetAtomicNum()
        if atom_num == 8:  # Oxygen
            if atom.GetTotalDegree() == 1:
                return True, "Molecule contains alcohol functional group"
            elif atom.GetTotalDegree() == 2:
                if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()):
                    return True, "Molecule contains aldehyde or ketone functional group"
                else:
                    return True, "Molecule contains ether functional group"
            elif atom.GetTotalDegree() == 3:
                return True, "Molecule contains carboxylic acid functional group"
        elif atom_num == 16:  # Sulfur
            return True, "Molecule contains sulfur functional group"
    
    return True, "Molecule matches the structural criteria for a monoterpenoid"