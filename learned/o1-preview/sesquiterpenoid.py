"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:26873 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is a terpenoid derived from a sesquiterpene (C15 hydrocarbons built from three isoprene units).
    The skeleton may be rearranged or modified by the removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12 or c_count > 15:
        return False, f"Carbon count is {c_count}, which is not in the range of 12-15"

    # Check for terpenoid functional groups (oxygen-containing groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms detected, not a terpenoid"

    # Estimate number of isoprene units by counting C=C bonds
    num_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_double_bonds +=1
    if num_double_bonds < 2:
        return False, f"Only {num_double_bonds} carbon-carbon double bonds detected, less than expected for sesquiterpenoids"

    # Check for isoprene unit patterns (heuristic)
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C=C")
    matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(matches) < 1:
        return False, "No isoprene units detected"

    # Evaluate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180 or mol_wt > 300:
        return False, f"Molecular weight of {mol_wt:.2f} is not typical for sesquiterpenoids"

    return True, "Molecule matches criteria for a sesquiterpenoid"