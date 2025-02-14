"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define SMARTS patterns for known biological molecules

    # Amino acid residue (general pattern)
    amino_acid_smarts = "[NX3;H2,H1;!$(N=C)][C;!$(C=O)](C(=O)[O,N])"
    amino_acid_pattern = Chem.MolFromSmarts(amino_acid_smarts)
    has_amino_acid = mol.HasSubstructMatch(amino_acid_pattern)
    if has_amino_acid:
        types_found.add('amino acid')

    # Peptide bond (O=C-N)
    peptide_bond_smarts = "C(=O)N"
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    has_peptide_bond = mol.HasSubstructMatch(peptide_bond_pattern)
    if has_peptide_bond:
        types_found.add('peptide')

    # Sugar (monosaccharide with multiple hydroxyl groups)
    sugar_smarts = "[C;R][O;R][C;R][C;R][C;R][C;R]"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    if has_sugar:
        types_found.add('sugar')

    # Nucleotide base (purine or pyrimidine rings)
    purine_smarts = "n1c[nH]c2c1ncnc2"
    pyrimidine_smarts = "c1cnc[nH]c1"
    purine_pattern = Chem.MolFromSmarts(purine_smarts)
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    if has_purine or has_pyrimidine:
        types_found.add('nucleotide')

    # Fatty acid (long aliphatic chain with terminal carboxylic acid)
    fatty_acid_smarts = "C(=O)[O;H1][C;X4][C;X4][C;X4][C;X4][C;X4]"
    fatty_acid_pattern = Chem.MolFromSmarts(fatty_acid_smarts)
    has_fatty_acid = mol.HasSubstructMatch(fatty_acid_pattern)
    if has_fatty_acid:
        types_found.add('fatty acid')

    # Coenzyme A (CoA) fragment (simplified pattern)
    coenzymeA_smarts = "C1=NC2=C(N1)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OCCNC(=O)CCNC(=O)C(O)C(C)(C)OP(=O)(O)O)O)O"
    coenzymeA_pattern = Chem.MolFromSmarts(coenzymeA_smarts)
    has_coenzymeA = mol.HasSubstructMatch(coenzymeA_pattern)
    if has_coenzymeA:
        types_found.add('coenzyme A')

    # Glutathione (tripeptide of glutamate, cysteine, and glycine)
    glutathione_smarts = "N[C@@H](CC[C@H](C(=O)O)N)C(=O)N[C@@H](CS)C(=O)NCC(=O)O"
    glutathione_pattern = Chem.MolFromSmarts(glutathione_smarts)
    has_glutathione = mol.HasSubstructMatch(glutathione_pattern)
    if has_glutathione:
        types_found.add('glutathione')

    # Determine if the molecule is a bioconjugate
    if len(types_found) >= 2:
        return True, f"Contains {', '.join(types_found)} covalently linked"
    else:
        return False, "Does not contain at least 2 different covalently linked biological units"