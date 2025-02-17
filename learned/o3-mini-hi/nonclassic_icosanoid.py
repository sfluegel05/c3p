"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: Nonclassic icosanoid
Definition: Any biologically active signalling molecule made by oxygenation of C20 fatty acids
            other than the classic icosanoids (the leukotrienes and the prostanoids).
This function uses heuristic structural filters:
  - The molecule should have exactly 20 carbon atoms (i.e. derived from a C20 fatty acid).
  - The molecule should be oxygenated (have several oxygen atoms).
  - The molecule should have several unsaturations (double bonds) indicative of a polyunsaturated
    fatty acid.
  - To exclude prostanoids, we check for a cyclopentane ring (a 5-membered non‐aromatic ring).
  - (Exclusion of leukotrienes by structure is challenging, so we rely on the above filters.)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    
    A nonclassic icosanoid is defined as an oxygenated derivative of a C20 fatty acid,
    excluding the classic icosanoids (leukotrienes and prostanoids). In our heuristic, the
    molecule must have 20 carbons, exhibit oxygenation and unsaturation typical for such a fatty
    acid derivative, and not contain a cyclopentane ring (a hallmark of prostanoids).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is judged to be a nonclassic icosanoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total carbons (atomic number 6)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Molecule has {c_count} carbon atoms; expected exactly 20 for a C20 fatty acid derivative."

    # Count total oxygen atoms (atomic number 8)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Molecule has insufficient oxygenation (only {o_count} oxygen atoms) to be an oxygenated fatty acid derivative."

    # Count double bonds (unsaturations)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 3:
        return False, f"Only {double_bond_count} double bonds found; typical icosanoids are polyunsaturated."

    # Check for the presence of a 5-membered ring 
    # (prostanoids feature a cyclopentane ring; excluding them is one filter for nonclassic molecules)
    ring_info = mol.GetRingInfo().AtomRings()
    for ring in ring_info:
        if len(ring) == 5:
            # Further check: if the ring is entirely non-aromatic (most cyclopentane rings are non-aromatic)
            # we assume that this ring may indicate a prostanoid.
            aromatic_atoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            if aromatic_atoms == 0:
                return False, "Molecule contains a cyclopentane ring, characteristic of prostanoids."

    # (Optional additional check: many icosanoids are acids or esters,
    # but nonclassic ones can also be methyl esters; we do not enforce this here.)

    # If all heuristic filters pass, then classify as a nonclassic icosanoid.
    return True, "Molecule is a C20 oxygenated fatty acid derivative lacking classic prostanoid/leukotriene features."