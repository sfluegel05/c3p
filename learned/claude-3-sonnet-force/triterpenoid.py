"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit import RDConfig
import os

# Load common triterpenoid skeletons
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseNameToMol.fsm')
suppl = Chem.SmilesMolSupplier(fdefName, titleLine=False)
skeletal_mols = [x for x in suppl if x is not None]

# Define SMARTS patterns for common modifications
REMOVE_METHYL = "[CH3]"
REARRANGE_RINGS = "[R2]@[r3,r4,r5,r6]@[R1]"
ADD_OXYGEN = "[OH],[O]"

def is_triterpenoid(smiles):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28 or c_count > 32:
        return False, f"Carbon count ({c_count}) outside the expected range for triterpenoids (28-32)"

    # Check for common triterpenoid skeletons
    mcs = rdFMCS.FindMCS(skeletal_mols + [mol], matchValences=True, completeRingsOnly=True, timeout=1)
    if mcs.numBonds == 0 or mcs.numAtoms == 0:
        return False, "No common triterpenoid skeleton found"

    # Check for common modifications
    remove_methyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(REMOVE_METHYL))
    rearrange_rings_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(REARRANGE_RINGS))
    add_oxygen_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(ADD_OXYGEN))

    if remove_methyl_matches or rearrange_rings_matches or add_oxygen_matches:
        return True, "Contains a triterpenoid skeleton with common modifications"

    # If no modifications found, check if it's a pristine triterpenoid skeleton
    if mcs.numAtoms == mol.GetNumAtoms() and mcs.numBonds == mol.GetNumBonds():
        return True, "Pristine triterpenoid skeleton"

    return False, "No evidence of being a triterpenoid"