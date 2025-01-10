"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid has a beta-D-glucuronic acid moiety connected via a glycosidic bond.

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

    # Define glucuronic acid pattern (beta-D-glucuronic acid with open glycosidic bond)
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)C(=O)O)O)O)O)[OX2]")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Verify the glycosidic bond is beta (C1-O-C')
    glycosidic_bond = None
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6 and
                begin_atom.GetDegree() == 2 and end_atom.GetDegree() >= 2):
                glycosidic_bond = bond
                break

    if glycosidic_bond is None:
        return False, "No glycosidic bond found"

    # Check if the glycosidic bond is beta (C1-O-C')
    glucuronic_acid_matches = mol.GetSubstructMatches(glucuronic_acid_pattern)
    for match in glucuronic_acid_matches:
        anomeric_carbon = match[0]
        oxygen = match[5]
        if (mol.GetBondBetweenAtoms(anomeric_carbon, oxygen) == glycosidic_bond and
            mol.GetAtomWithIdx(anomeric_carbon).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            return True, "Contains beta-D-glucuronic acid moiety connected via glycosidic bond"

    return False, "No beta-glycosidic bond found"