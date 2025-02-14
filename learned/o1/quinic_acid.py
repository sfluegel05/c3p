"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or a derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid - a cyclohexane ring with multiple hydroxyl groups
    (or their derivatives like esters) and a carboxylic acid or ester group attached to the ring.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is quinic acid or derivative, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in the molecule"

    # Define patterns
    carboxylic_acid_smarts = "[C](=O)[O;H1,H2]"  # Carboxylic acid group
    ester_smarts = "[C](=O)[O;!H]"               # Ester group
    hydroxyl_smarts = "[OX2H]"                   # Hydroxyl group
    ether_smarts = "[OX2H0]"                     # Ether linkage (oxygen without hydrogen)

    carboxylic_acid = Chem.MolFromSmarts(carboxylic_acid_smarts)
    ester = Chem.MolFromSmarts(ester_smarts)
    hydroxyl = Chem.MolFromSmarts(hydroxyl_smarts)
    ether = Chem.MolFromSmarts(ether_smarts)

    # Iterate over all six-membered rings
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Skip non-six-membered rings

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check if all atoms in the ring are sp3 carbons (saturated cyclohexane)
        if not all(atom.GetAtomicNum() == 6 and atom.GetDegree() == 4 for atom in ring_atoms):
            continue  # Skip if ring is not a cyclohexane

        substituent_count = 0
        carboxylic_acid_found = False

        # Check substituents on ring atoms
        for atom in ring_atoms:
            atom_idx = atom.GetIdx()
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring:
                    continue  # Skip ring atoms

                # Check for carboxylic acid or ester group
                bond_atoms = {atom_idx, nbr_idx}
                bond_smarts = Chem.MolFragmentToSmiles(mol, atomsToUse=list(bond_atoms), canonical=False)
                sub_mol = Chem.MolFromSmiles(bond_smarts)

                if sub_mol.HasSubstructMatch(carboxylic_acid):
                    carboxylic_acid_found = True
                elif sub_mol.HasSubstructMatch(ester):
                    carboxylic_acid_found = True  # Considering esters of quinic acid
                elif nbr.GetAtomicNum() == 8:
                    # Check for hydroxyl or ether groups
                    if nbr.GetDegree() == 1:
                        # Hydroxyl group
                        substituent_count += 1
                    elif nbr.GetDegree() == 2:
                        # Possible ether or ester linkage
                        substituent_count += 1

        if substituent_count >= 3 and carboxylic_acid_found:
            return True, "Contains cyclitol carboxylic acid core (quinic acid or derivative)"
        else:
            reasons = []
            if substituent_count < 3:
                reasons.append(f"Found {substituent_count} hydroxyl or ester groups on ring, need at least 3")
            if not carboxylic_acid_found:
                reasons.append("No carboxylic acid or ester group attached to ring")
            return False, "; ".join(reasons)

    return False, "No suitable cyclohexane ring found with required substituents"