"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion.
    The structure must have:
    - One carboxylic acid group (deprotonated as -COO-)
    - Oxo group (=O) at the alpha (2nd) position
    - Overall negative charge
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for exactly one carboxylic acid group (as anion)
    carboxylate = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups (needs 1)"

    # Get the carboxylic acid carbon
    carboxyl_carbon = carboxylate_matches[0][0]

    # Check adjacent carbon (alpha position) has an oxo group
    alpha_carbon = None
    for bond in mol.GetAtomWithIdx(carboxyl_carbon).GetBonds():
        neighbor = bond.GetOtherAtomIdx(carboxyl_carbon)
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6:  # Carbon
            alpha_carbon = neighbor
            break

    if alpha_carbon is None:
        return False, "No alpha carbon found"

    # Check alpha carbon has an oxo group (=O)
    alpha_oxo = False
    for bond in mol.GetAtomWithIdx(alpha_carbon).GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            neighbor = bond.GetOtherAtomIdx(alpha_carbon)
            if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 8:  # Oxygen
                alpha_oxo = True
                break

    if not alpha_oxo:
        return False, "Alpha carbon lacks oxo group"

    # Check overall charge is negative (anion)
    charge = Chem.GetFormalCharge(mol)
    if charge >= 0:
        return False, f"Charge is {charge}, must be negative"

    return True, "2-oxo monocarboxylic acid anion"