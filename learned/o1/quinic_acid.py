"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: quinic acid
"""
from rdkit import Chem

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
    carboxylic_acid_or_ester_smarts = "[CX3](=O)[OX1H0,R0]"
    oxygen_substituent_smarts = "[#6][OX2H0,H1]"  # Carbon attached to oxygen (hydroxyl, ether, ester)

    carboxylic_acid_or_ester = Chem.MolFromSmarts(carboxylic_acid_or_ester_smarts)
    oxygen_substituent = Chem.MolFromSmarts(oxygen_substituent_smarts)

    # Iterate over all six-membered rings
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Skip non-six-membered rings

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        substituent_oxygens = set()
        carboxylic_acid_or_ester_found = False

        # Check substituents on ring atoms
        for atom in ring_atoms:
            # Check for oxygen substituents attached to ring carbons
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx in ring:
                    continue  # Skip atoms within the ring

                # Check if the neighbor is an oxygen atom
                if neighbor.GetAtomicNum() == 8:
                    substituent_oxygens.add(nbr_idx)
                else:
                    # Check for carboxylic acid or ester group attached to the ring carbon
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue  # We are interested in single bonds only
                    substructure = Chem.PathToSubmol(mol, [atom.GetIdx(), neighbor.GetIdx()])
                    if substructure.HasSubstructMatch(carboxylic_acid_or_ester):
                        carboxylic_acid_or_ester_found = True

        # Check if ring meets the criteria
        if len(substituent_oxygens) >= 3 and carboxylic_acid_or_ester_found:
            return True, "Contains cyclitol carboxylic acid core (quinic acid or derivative)"
        else:
            reasons = []
            if len(substituent_oxygens) < 3:
                reasons.append(f"Found {len(substituent_oxygens)} oxygen substituents on ring carbons, need at least 3")
            if not carboxylic_acid_or_ester_found:
                reasons.append("No carboxylic acid or ester group attached to ring carbons")
            return False, "; ".join(reasons)

    return False, "No suitable cyclohexane ring found with required substituents"