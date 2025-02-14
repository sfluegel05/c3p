"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is defined as a glucosiduronic acid resulting from the
    formal condensation of any substance with beta-D-glucuronic acid to form a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define beta-D-glucuronic acid substructure
    beta_d_glucuronic_acid = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)C(=O)O")

    # Check for beta-D-glucuronic acid substructure
    matches = mol.GetSubstructMatches(beta_d_glucuronic_acid)
    if not matches:
        return False, "No beta-D-glucuronic acid substructure found"

    # Define glycosidic bond pattern
    glycosidic_bond = Chem.MolFromSmarts("[OX2]C")

    # Check for glycosidic bond connecting beta-D-glucuronic acid to the rest of the molecule
    for match in matches:
        o_idx = match[0]  # Index of the ring oxygen in beta-D-glucuronic acid
        neighbors = [mol.GetBondBetweenAtoms(o_idx, n).GetBondType() for n in mol.GetAtomWithIdx(o_idx).GetNeighbors()]
        if Chem.BondType.SINGLE in neighbors and mol.HasSubstructMatch(glycosidic_bond, atomIds=[o_idx]):
            # Additional checks (optional)
            mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 300:  # Adjust as needed
                return False, "Molecular weight too low for beta-D-glucosiduronic acid"

            return True, "Contains beta-D-glucuronic acid substructure connected via a glycosidic bond"

    return False, "No glycosidic bond connecting beta-D-glucuronic acid found"