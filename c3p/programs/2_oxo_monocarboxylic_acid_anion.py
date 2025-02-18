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
    - Exactly one carboxylic acid group (deprotonated as -COO-)
    - Oxo group (=O) at the alpha (2nd) position as a ketone (not amide, ester, etc.)
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

    carboxyl_carbon = carboxylate_matches[0][0]

    # Find alpha carbon (adjacent to carboxyl carbon)
    alpha_carbons = []
    for bond in mol.GetAtomWithIdx(carboxyl_carbon).GetBonds():
        neighbor = bond.GetOtherAtomIdx(carboxyl_carbon)
        if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6:  # Carbon
            alpha_carbons.append(neighbor)

    if not alpha_carbons:
        return False, "No alpha carbon found"

    # Check each alpha carbon for ketone group
    valid_alpha = False
    for alpha in alpha_carbons:
        alpha_atom = mol.GetAtomWithIdx(alpha)
        # Check for double bond to oxygen
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtomIdx(alpha)
                neighbor_atom = mol.GetAtomWithIdx(neighbor)
                if neighbor_atom.GetAtomicNum() == 8:  # Oxygen
                    # Check oxygen has no other bonds
                    if neighbor_atom.GetDegree() == 1:
                        # Check other neighbors of alpha are carbons
                        valid = True
                        for b in alpha_atom.GetBonds():
                            other_idx = b.GetOtherAtomIdx(alpha)
                            if other_idx == neighbor:  # skip the oxygen
                                continue
                            other_atom = mol.GetAtomWithIdx(other_idx)
                            if other_atom.GetAtomicNum() != 6:
                                valid = False
                                break
                        if valid:
                            valid_alpha = True
                            break
        if valid_alpha:
            break

    if not valid_alpha:
        return False, "Alpha carbon lacks ketone group or has non-carbon substituents"

    # Check overall charge is negative
    charge = Chem.GetFormalCharge(mol)
    if charge >= 0:
        return False, f"Charge is {charge}, must be negative"

    return True, "2-oxo monocarboxylic acid anion"