"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove hydrogens for simplicity
    mol = Chem.RemoveHs(mol)

    # Check for triterpenoid core (C30 skeleton)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 30:
        return False, f"Contains only {num_carbons} carbon atoms; fewer than 30 needed for triterpenoid"

    # Identify ring systems
    ssr = Chem.GetSSSR(mol)
    if ssr < 4:
        return False, f"Contains only {ssr} rings; triterpenoids typically have multiple rings"
    
    # Find the largest ring system (assumed to be the core)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    largest_ring_system = max(atom_rings, key=len)
    core_atoms = set(largest_ring_system)

    # Detect sugar moieties
    # Define SMARTS patterns for sugars (pyranose and furanose rings)
    sugar_patterns = [
        '[OX2H][CX4H]',  # Hydroxyl group attached to carbon
        '[CX4H1][OX2H][CX4]'  # Ether linkage pattern
    ]
    sugar_mols = [Chem.MolFromSmarts(pat) for pat in sugar_patterns]
    if not all(sugar_mols):
        return False, "Failed to parse sugar SMARTS patterns"

    # Find sugar atoms
    sugar_atoms = set()
    for sugar_mol in sugar_mols:
        matches = mol.GetSubstructMatches(sugar_mol)
        for match in matches:
            sugar_atoms.update(match)

    if not sugar_atoms:
        return False, "No sugar moieties found"

    # Check for glycosidic linkage between core and sugar
    found_linkage = False
    for bond in mol.GetBonds():
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()

        # Check for C-O-C linkage between core and sugar
        if bond.GetBondType() == Chem.BondType.SINGLE:
            if (atom1_idx in core_atoms and atom2_idx in sugar_atoms) or \
               (atom2_idx in core_atoms and atom1_idx in sugar_atoms):
                if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8) or \
                   (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6):
                    found_linkage = True
                    break

    if not found_linkage:
        return False, "No glycosidic linkage between core and sugar found"

    # Optional: Check for common functional groups in triterpenoid saponins
    hydroxyls = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1]
    carboxyls = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyls)

    if len(hydroxyls) < 1 and len(carboxyl_matches) < 1:
        return False, "No hydroxyl or carboxyl groups found; unusual for triterpenoid saponins"

    return True, "Molecule is a triterpenoid saponin with triterpenoid core and glycosidically linked sugar moiety"