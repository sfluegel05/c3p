"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is a molecular entity consisting of at least 2 different biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    types_found = set()

    # Define mol objects for known biological molecules

    # Glutathione
    glutathione_smiles = "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"
    glutathione_mol = Chem.MolFromSmiles(glutathione_smiles)
    has_glutathione = mol.HasSubstructMatch(glutathione_mol)
    if has_glutathione:
        types_found.add('glutathione')

    # Coenzyme A (simplified pattern)
    coenzymeA_smiles = "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO)C(O)C(O)C1O)n2cnc3c(N)ncnc23"
    coenzymeA_mol = Chem.MolFromSmiles(coenzymeA_smiles)
    has_coenzymeA = mol.HasSubstructMatch(coenzymeA_mol)
    if has_coenzymeA:
        types_found.add('coenzyme A')

    # Fatty acid chain (simplified long chain carboxylic acid)
    fatty_acid_smarts = "C(=O)[CX4;!$(*=,#[!#6])]CCCC"  # Carboxyl group attached to a long chain
    fatty_acid_pattern = Chem.MolFromSmarts(fatty_acid_smarts)
    has_fatty_acid = mol.HasSubstructMatch(fatty_acid_pattern)
    if has_fatty_acid:
        types_found.add('fatty acid')

    # Sugar (simplified glucose ring)
    sugar_smarts = "C1(CO)OC(O)C(O)C(O)C1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    if has_sugar:
        types_found.add('sugar')

    # Nucleotide (simplified purine/pyrimidine ring with sugar and phosphate)
    nucleotide_smarts = "n1c([nH])cnc1C2OC(COP(O)(O)=O)C(O)C2O"
    nucleotide_pattern = Chem.MolFromSmarts(nucleotide_smarts)
    has_nucleotide = mol.HasSubstructMatch(nucleotide_pattern)
    if has_nucleotide:
        types_found.add('nucleotide')

    # Amino acid residue (exclude peptide bonds to prevent overcounting)
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")  # General amino acid pattern
    has_amino_acid = mol.HasSubstructMatch(amino_acid_pattern)
    if has_amino_acid:
        types_found.add('amino acid')

    # Peptide bond pattern to detect peptides
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)NC")  # Peptide bond
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # If only amino acids are found and the molecule has peptide bonds, consider it a peptide (single unit)
    if types_found == {'amino acid'} and num_peptide_bonds >= 1:
        types_found = {'peptide'}

    # Determine if the molecule is a bioconjugate
    if len(types_found) >= 2:
        return True, f"Contains {', '.join(types_found)} covalently linked"
    else:
        return False, "Does not contain at least 2 different covalently linked biological units"