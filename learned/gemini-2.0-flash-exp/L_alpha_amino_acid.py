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

    # SMARTS pattern to match L-alpha-amino acid core structure with L-configuration
    #   [#6] : Carbon atom
    #   [C@@H] : Carbon atom with L-stereoconfiguration
    #   [NX3;H2,H3;!+] : Nitrogen with three bonds (either -NH2 or -NH3+)
    #   [CX3](=[OX1]) : Carboxyl group (-COOH or -COO-)
    #   ~ : Indicates any kind of bond
    #   [!#6] any atom that is not C (for the carboxyl)
    #   [!#1] any atom that is not hydrogen (for the amino group)
    #   [!$(C(=O)N)] the carboxyl group should not be bonded to a nitrogen atom (prevents matching peptides)
    #   [C;!$(C(=O)N)] The carbon attached to the alpha carbon can not be part of a peptide bond.
    #   [!$(NC(=O))] The N should not be part of a peptide bond
    #   [C@H1] match exactly one hydrogen

    # SMARTS for the L-alpha-amino acid core (neutral and zwitterionic forms, with more restrictions)
    smarts_neutral = "[NX3;H2;!+][C@@H1]([!#1;!$(C(=O)N)][CX3](=[OX1])[!#6;!$(NC(=O))])"
    smarts_zwitter = "[NX3;H3;!+][C@@H1]([!#1;!$(C(=O)N)][CX3](=[OX1])[OX1-])"
    smarts_zwitter2 = "[NX3;H3;!+][C@@H1]([!#1;!$(C(=O)N)][CX3](=O)[O-])"


    # SMARTS for cyclic L-alpha-amino acids
    smarts_cyclic_neutral =  "[NX3;H1;!+][C@@H1]([C;!$(C(=O)N)]-[CX3](=[OX1])[!#6;!$(NC(=O))])" # alpha C is within a ring
    smarts_cyclic_zwitter =  "[NX3;H2;!+][C@@H1]([C;!$(C(=O)N)]-[CX3](=[OX1])[OX1-])"
    smarts_cyclic_zwitter2 = "[NX3;H2;!+][C@@H1]([C;!$(C(=O)N)]-[CX3](=O)[O-])"


    # Check if the molecule matches the SMARTS pattern
    core_pattern_neutral = Chem.MolFromSmarts(smarts_neutral)
    core_pattern_zwitter = Chem.MolFromSmarts(smarts_zwitter)
    core_pattern_zwitter2 = Chem.MolFromSmarts(smarts_zwitter2)
    core_pattern_cyclic_neutral = Chem.MolFromSmarts(smarts_cyclic_neutral)
    core_pattern_cyclic_zwitter = Chem.MolFromSmarts(smarts_cyclic_zwitter)
    core_pattern_cyclic_zwitter2 = Chem.MolFromSmarts(smarts_cyclic_zwitter2)


    if not (mol.HasSubstructMatch(core_pattern_neutral) or
             mol.HasSubstructMatch(core_pattern_zwitter) or
             mol.HasSubstructMatch(core_pattern_zwitter2) or
             mol.HasSubstructMatch(core_pattern_cyclic_neutral) or
             mol.HasSubstructMatch(core_pattern_cyclic_zwitter) or
             mol.HasSubstructMatch(core_pattern_cyclic_zwitter2)
             ):
         return False, "Molecule does not contain the L-alpha-amino acid core structure."

    # count the number of carbonyl groups. Must be exactly one to be a single amino acid and not a peptide or a molecule with additional carbonyls
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) != 1:
        return False, f"Molecule contains {len(carbonyl_matches)} carbonyl groups. Must contain exactly one."

    # if it matches the L-alpha-amino acid structure, we return True
    return True, "Molecule contains the L-alpha-amino acid core structure with L-configuration."