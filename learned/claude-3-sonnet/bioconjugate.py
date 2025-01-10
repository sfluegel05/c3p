"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
A molecular entity consisting of at least 2 biological molecules covalently linked together.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate must contain at least 2 biological molecules covalently linked.

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

    # Initialize counters for biological components
    bio_components = 0
    reasons = []

    # Check for amino acid patterns
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        amino_acid_matches = len(mol.GetSubstructMatches(amino_acid_pattern))
        if amino_acid_matches >= 1:
            bio_components += 1
            reasons.append(f"Contains {amino_acid_matches} amino acid(s)")

    # Check for peptide bonds
    peptide_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(peptide_pattern):
        bio_components += 1
        reasons.append("Contains peptide bond")

    # Check for nucleotide/CoA patterns
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")  # Adenine base
    if mol.HasSubstructMatch(nucleotide_pattern):
        bio_components += 1
        reasons.append("Contains nucleotide base")

    # Check for sugar patterns
    sugar_pattern = Chem.MolFromSmarts("[OH1][CX4H]1[OX2][CX4H]([CX4H]([OH1])[CX4H]([OH1])[CX4H]1[OH1])")
    if mol.HasSubstructMatch(sugar_pattern):
        bio_components += 1
        reasons.append("Contains sugar moiety")

    # Check for fatty acid patterns (long carbon chain with terminal acid)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]C(=O)[OH]")
    if mol.HasSubstructMatch(fatty_acid_pattern):
        bio_components += 1
        reasons.append("Contains fatty acid chain")

    # Check for glutathione pattern
    glutathione_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H](CCC(=O)[NX3,NX4+][CX4H](CS)[CX3](=O)[NX3,NX4+][CX4H])[CX3](=O)[O-,OH]")
    if mol.HasSubstructMatch(glutathione_pattern):
        bio_components += 1
        reasons.append("Contains glutathione")

    # Check for CoA pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[CH]([OH])C(C)(C)COP([OH])(=O)OP([OH])(=O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP([OH])([OH])=O")
    if mol.HasSubstructMatch(coa_pattern):
        bio_components += 1
        reasons.append("Contains Coenzyme A")

    # Check for thioether linkage (common in conjugates)
    thioether_pattern = Chem.MolFromSmarts("[CX4]S[CX4]")
    if mol.HasSubstructMatch(thioether_pattern):
        reasons.append("Contains thioether linkage")

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for bioconjugate"

    # Final decision
    if bio_components >= 2:
        return True, "Bioconjugate: " + "; ".join(reasons)
    elif bio_components == 1:
        return False, "Only one biological component found: " + "; ".join(reasons)
    else:
        return False, "No biological components identified"