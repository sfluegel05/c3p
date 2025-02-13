"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:36035 tetraterpenoid
A tetraterpenoid is any terpenoid derived from a tetraterpene. The term includes compounds in which the C40
skeleton of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms
(generally methyl groups).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdqueries
from typing import Tuple

def is_tetraterpenoid(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check number of carbon atoms (typically 40-60)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 40 or n_carbons > 60:
        return False, f"Number of carbon atoms ({n_carbons}) outside expected range for tetraterpenoids (40-60)"

    # Check molecular weight (typically 500-1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1000:
        return False, f"Molecular weight ({mol_wt:.2f} Da) outside expected range for tetraterpenoids (500-1000 Da)"

    # Check for presence of isoprene units
    isoprene_pattern = Chem.MolFromSmarts("[C@H]([C@H](C)[C@H](C)C)=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 8:
        return False, "Fewer than 8 isoprene units detected"

    # Check for presence of long carbon chains and/or rings
    chain_pattern = Chem.MolFromSmarts("CCC(C)CCC")
    ring_pattern = rdqueries.IsRingQuery()
    has_long_chains = bool(mol.GetSubstructMatches(chain_pattern))
    has_rings = any(ring_pattern(ring) for ring in mol.GetRingInfo().AtomRings())
    if not has_long_chains and not has_rings:
        return False, "No long carbon chains or rings detected"

    return True, "Structure consistent with a tetraterpenoid (derived from a tetraterpene with C40 skeleton)"