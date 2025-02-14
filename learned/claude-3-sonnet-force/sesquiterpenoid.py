"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:36689 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is a terpenoid derived from a sesquiterpene skeleton (C15),
    which may be rearranged or modified by the removal of atoms.

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

    # Check molecular weight range (typically 200-300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 300:
        return False, f"Molecular weight {mol_wt:.2f} outside typical range for sesquiterpenoids"

    # Count atoms and check for C15 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Less than 15 carbon atoms, not a sesquiterpene skeleton"

    # Check for common functional groups (alcohols, carbonyls, etc.)
    has_alcohol = any(atom.GetHybridization() == Chem.HybridizationAtomType.SP3 and
                      atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    has_carbonyl = any(atom.GetHybridization() == Chem.HybridizationAtomType.SP2 and
                       atom.GetIsAromatic() and atom.GetFormalCharge() == 0 for atom in mol.GetAtoms())

    if not has_alcohol and not has_carbonyl:
        return False, "No common functional groups found for sesquiterpenoids"

    # Check for rings and unsaturated bonds
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_unsat_bonds = rdMolDescriptors.CalcNumUnsaturatedCarbonAtoms(mol)

    if n_rings == 0 and n_unsat_bonds == 0:
        return False, "No rings or unsaturated bonds, unlikely to be a sesquiterpenoid"

    return True, "Meets criteria for a sesquiterpenoid: C15 skeleton, functional groups, rings/unsaturation"