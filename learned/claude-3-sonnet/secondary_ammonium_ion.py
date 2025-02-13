"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:51953 secondary ammonium ion
An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if any nitrogen atom is a secondary ammonium ion
    for nitrogen in nitrogen_atoms:
        # Check if the nitrogen has positive charge and sp3 hybridization
        if nitrogen.GetFormalCharge() == 1 and nitrogen.GetHybridization() == Chem.HybridizationType.SP3:
            # Check if the nitrogen has exactly two substituents (excluding hydrogens)
            non_hydrogen_neighbors = [neighbor for neighbor in nitrogen.GetNeighbors() if neighbor.GetAtomicNum() != 1]
            if len(non_hydrogen_neighbors) == 2:
                # Check if the nitrogen is part of an amine group
                neighboring_groups = [get_neighboring_group(mol, nitrogen, neighbor) for neighbor in non_hydrogen_neighbors]
                if any(group.GetAtomicNum() == 6 for group in neighboring_groups):
                    return True, "Molecule contains a secondary ammonium ion group"

    return False, "No secondary ammonium ion group found"

def get_neighboring_group(mol, atom, neighbor):
    """
    Finds the group connected to the given atom through the specified neighbor.

    Args:
        mol (Mol): RDKit Mol object
        atom (Atom): RDKit Atom object
        neighbor (Atom): RDKit Atom object representing the neighbor connecting the group

    Returns:
        AtomPtrContainer: RDKit container of atoms in the neighboring group
    """
    env = Chem.FindAtomEnvironmentOfRadiusN(mol, neighbor, 2)
    atoms = env.GetAtomPtrs()
    group = Chem.AtomPtrContainer()
    for atom in atoms:
        if atom != neighbor and atom != atom:
            group.Add(atom)
    return group