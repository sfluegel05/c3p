"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35701 dicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate tautomers and resonance structures
    tautomers = list(AllChem.ResonanceMolSupplier(mol, AllChem.KIDERS_ALL))
    tautomers.append(mol)

    # Check for two carboxyl groups
    carboxyl_pattern = Chem.MolFromSmarts("[C](=[O])(O)")
    for tautomer in tautomers:
        carboxyl_matches = tautomer.GetSubstructMatches(carboxyl_pattern)
        if len(carboxyl_matches) >= 2:
            # Check if the carboxyl groups are part of larger functional groups like esters or amides
            ester_pattern = Chem.MolFromSmarts("[C](=[O])(O)[O]")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])(O)[N]")
            if not any(tautomer.HasSubstructMatch(ester_pattern) or tautomer.HasSubstructMatch(amide_pattern)):
                return True, "Contains two or more carboxyl groups not part of esters or amides"

    # Check for specific substructures or functional groups commonly found in dicarboxylic acids
    alpha_keto_acid_pattern = Chem.MolFromSmarts("[C](=[O])(O)[C](=[O])")
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2]([C]([C](=[O])(O))([C](=[O])(O)))")
    cyclic_dicarboxylic_acid_pattern = Chem.MolFromSmarts("[C]([C](=[O])(O))([C](=[O])(O))1[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]1([H])([H])")

    for tautomer in tautomers:
        if tautomer.HasSubstructMatch(alpha_keto_acid_pattern) or \
           tautomer.HasSubstructMatch(amino_acid_pattern) or \
           tautomer.HasSubstructMatch(cyclic_dicarboxylic_acid_pattern):
            return True, "Contains specific substructures or functional groups commonly found in dicarboxylic acids"

    return False, "Does not contain two or more carboxyl groups or specific dicarboxylic acid substructures"