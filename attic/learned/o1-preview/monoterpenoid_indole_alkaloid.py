"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety derived from L-tryptophan and
    a monoterpene-derived moiety (C10 unit), usually linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for indole moiety
    indole_pattern = Chem.MolFromSmarts('c1c([n&H0])ccc1')  # Five-membered ring with nitrogen
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')  # Six-membered benzene ring

    # Check for fused indole (benzene ring fused to pyrrole ring)
    fused_indole = False
    ri = mol.GetRingInfo()

    # Look for fused rings
    for ring1 in ri.AtomRings():
        if len(ring1) == 5:
            atoms_in_ring1 = set(ring1)
            # Check if ring contains nitrogen
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring1):
                for ring2 in ri.AtomRings():
                    if len(ring2) == 6:
                        atoms_in_ring2 = set(ring2)
                        # Check if rings are fused (share two atoms)
                        if len(atoms_in_ring1 & atoms_in_ring2) >= 2:
                            fused_indole = True
                            break
        if fused_indole:
            break

    if not fused_indole:
        return False, "No indole moiety found (fused benzene and pyrrole rings)"

    # Check for monoterpene-derived moiety (approximate)
    # Monoterpenes are C10 units built from isoprene units (C5H8)
    # Look for isoprene units or count number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:  # Adjusted to account for the indole ring and monoterpene unit
        return False, f"Too few carbons ({num_carbons}) to be a monoterpenoid indole alkaloid"

    # Optionally, check for alkaloid nitrogen (basic nitrogen atom)
    basic_nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 0 and atom.GetHybridization() != Chem.HybridizationType.SP2]
    if not basic_nitrogens:
        return False, "No basic nitrogen atom found typical of alkaloids"

    return True, "Contains fused indole moiety and sufficient size to be a monoterpenoid indole alkaloid"

__metadata__ = {
    'chemical_class': {
        'name': 'monoterpenoid indole alkaloid',
        'definition': 'A terpenoid indole alkaloid which is biosynthesised from L-tryptophan and diisoprenoid (usually secologanin) building blocks.'
    },
    'config': {
        # Configuration parameters if any
    }
}