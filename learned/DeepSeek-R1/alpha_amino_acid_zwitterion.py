"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:59999 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion has a protonated amino group ([NH3+]) and a deprotonated carboxyl group ([O-]),
    with both groups attached to adjacent carbons (alpha carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find ammonium group ([NH3+]) attached to a carbon
    ammonium_pattern = Chem.MolFromSmarts("[NH3+]")
    ammonium_matches = mol.GetSubstructMatches(ammonium_pattern)
    if not ammonium_matches:
        return False, "No ammonium group ([NH3+]) found"

    # Find carboxylate group ([O-]C(=O))
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group ([O-]C(=O)) found"

    # Check if ammonium and carboxylate are on adjacent carbons (alpha-amino structure)
    # Get the atoms involved in each group
    for ammonium_atom in ammonium_matches:
        for carboxylate_atom in carboxylate_matches:
            # Carboxylate atom is the carbon adjacent to the O-
            carboxyl_carbon = carboxylate_atom[1]
            # Ammonium atom is the nitrogen
            nitrogen = ammonium_atom[0]
            # Get the alpha carbon (attached to nitrogen)
            alpha_carbon = None
            for neighbor in mol.GetAtomWithIdx(nitrogen).GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    alpha_carbon = neighbor.GetIdx()
                    break
            if alpha_carbon is None:
                continue  # No adjacent carbon to nitrogen

            # Check if this alpha carbon is adjacent to the carboxyl carbon
            alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon)
            for neighbor in alpha_carbon_atom.GetNeighbors():
                if neighbor.GetIdx() == carboxyl_carbon:
                    return True, "Ammonium and carboxylate groups on adjacent carbons (alpha-amino zwitterion)"

    return False, "Ammonium and carboxylate groups not properly positioned for alpha-amino zwitterion"