"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.
    Quinic acid is characterized by a cyclohexane ring with multiple hydroxyl groups and a
    carboxylic acid group attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid or its derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to ensure proper valence perception
    try:
        Chem.Kekulize(mol)
    except:
        pass

    # Find cyclohexane rings (6-membered saturated rings)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    cyclohexane_rings = []
    for ring in rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetHybridization() == Chem.HybridizationType.SP3 for atom in atoms_in_ring):
                cyclohexane_rings.append(ring)

    if not cyclohexane_rings:
        return False, "No cyclohexane ring found"

    # For each cyclohexane ring, check for hydroxyl/ester groups and carboxylic acid
    for ring in cyclohexane_rings:
        ring_atoms = set(ring)
        hydroxyl_count = 0
        esterified_count = 0
        carboxylic_acid_found = False

        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check for oxygen-containing substituents (excluding ring oxygens)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx not in ring_atoms:
                    if neighbor.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            hydroxyl_count += 1  # Hydroxyl group
                    elif neighbor.GetAtomicNum() == 6:
                        # Check for ester group (O-C(=O)-)
                        esters = Chem.MolFromSmarts('[$(OC(=O)[#6])]')
                        if mol.HasSubstructMatch(esters):
                            esterified_count += 1

            # Check for carboxylic acid attached to the ring
            carboxylic_acid = Chem.MolFromSmarts('[C;R1][C](=O)[O;H1]')
            if mol.HasSubstructMatch(carboxylic_acid):
                carboxylic_acid_found = True

        total_oxygen_substituents = hydroxyl_count + esterified_count
        if total_oxygen_substituents < 3:
            continue  # Try the next ring
        if not carboxylic_acid_found:
            continue  # Try the next ring

        # Passed all checks
        return True, "Contains cyclohexane ring with multiple hydroxyl/ester groups and a carboxylic acid"

    return False, "Does not match quinic acid structural criteria"