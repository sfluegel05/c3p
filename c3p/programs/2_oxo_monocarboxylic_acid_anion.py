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
    - Exactly one carboxylic acid group (deprotonated as -COO⁻)
    - Oxo group (=O) at the alpha (2nd) position as a ketone
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

    # Find alpha carbon(s) (directly adjacent to carboxyl carbon)
    alpha_carbons = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors()
                     if neighbor.GetAtomicNum() == 6]

    if not alpha_carbons:
        return False, "No alpha carbon adjacent to carboxyl group"

    # Check if any alpha carbon has a ketone group
    ketone_found = False
    for alpha_idx in alpha_carbons:
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # Check for double bond to oxygen
        for bond in alpha_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtomIdx(alpha_idx)
                neighbor_atom = mol.GetAtomWithIdx(neighbor)
                if neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetDegree() == 1:
                    ketone_found = True
                    break
        if ketone_found:
            break

    if not ketone_found:
        return False, "No ketone group found on alpha carbon"

    # Verify overall negative charge
    charge = Chem.GetFormalCharge(mol)
    if charge >= 0:
        return False, f"Charge is {charge}, must be negative"

    return True, "2-oxo monocarboxylic acid anion"