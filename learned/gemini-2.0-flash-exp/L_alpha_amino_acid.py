"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    A single L-alpha-amino acid should have a single chiral center of the L configuration,
    a single amino group (-NH2 or -NH3+) and a single carboxyl group (-COOH or -COO-) directly linked to that chiral center.
    It should not be part of a larger peptide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the core L-alpha-amino acid structure.
    # Includes neutral and zwitterionic forms, as well as ring structures.
    # Specific constraints are applied for the alpha carbon, nitrogen, and carboxylic acid.
    # -  [NX3;H0,H1,H2,H3;!+]: Nitrogen with 0-3 explicit hydrogens, not charged, with any number of implicit H.
    # -  [C@@H1]: chiral carbon with explicit hydrogen and L stereochemistry
    # - [CX3](=[OX1])O: neutral carboxyl group (-COOH)
    # - [CX3](=O)[O-]: anionic carboxyl group (-COO-)
    # - [!#1]: any atom not being hydrogen. This ensures that the alpha carbon is only bonded to one non-H group other than NH and COOH

    smarts_neutral = "[NX3;H0,H1,H2,H3;!+][C@@H1]([CX4][CX3](=[OX1])O)[!#1]"
    smarts_zwitter = "[NX3;H0,H1,H2,H3;!+][C@@H1]([CX4][CX3](=O)[O-])[!#1]"

    smarts_cyclic_neutral =  "[NX3;H0,H1,H2,H3;!+][C@@H1;R]([CX4;R][CX3](=[OX1])O)[!#1]"
    smarts_cyclic_zwitter =  "[NX3;H0,H1,H2,H3;!+][C@@H1;R]([CX4;R][CX3](=O)[O-])[!#1]"
    
    
    # Check if molecule matches any of the core patterns
    core_pattern_neutral = Chem.MolFromSmarts(smarts_neutral)
    core_pattern_zwitter = Chem.MolFromSmarts(smarts_zwitter)
    core_pattern_cyclic_neutral = Chem.MolFromSmarts(smarts_cyclic_neutral)
    core_pattern_cyclic_zwitter = Chem.MolFromSmarts(smarts_cyclic_zwitter)
    
    if not (mol.HasSubstructMatch(core_pattern_neutral) or
             mol.HasSubstructMatch(core_pattern_zwitter) or
             mol.HasSubstructMatch(core_pattern_cyclic_neutral) or
             mol.HasSubstructMatch(core_pattern_cyclic_zwitter)):
         return False, "Molecule does not contain the L-alpha-amino acid core structure."


    # Check that there is only one carboxylic acid group present.
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX1H0-,OX2-]") # Matches -COOH or -COO-
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 1:
         return False, f"Molecule contains {len(carboxyl_matches)} carboxylic groups. Must contain exactly one."


    return True, "Molecule contains the L-alpha-amino acid core structure with L-configuration."