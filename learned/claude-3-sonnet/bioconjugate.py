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
    components = []

    # Check for peptide/amino acid components (count as one component)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H]([*])[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        amino_acid_matches = len(mol.GetSubstructMatches(amino_acid_pattern))
        if amino_acid_matches >= 1:
            components.append(f"peptide/amino acid chain ({amino_acid_matches} residues)")

    # Check for nucleotide components
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")  # Adenine base
    if mol.HasSubstructMatch(nucleotide_pattern):
        components.append("nucleotide")

    # Check for sugar components (excluding those part of nucleotides)
    sugar_pattern = Chem.MolFromSmarts("[OH1][CX4H]1[OX2][CX4H]([CX4H]([OH1])[CX4H]([OH1])[CX4H]1[OH1])")
    if mol.HasSubstructMatch(sugar_pattern):
        components.append("sugar")

    # Check for fatty acid components
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]C(=O)[OH]")
    if mol.HasSubstructMatch(fatty_acid_pattern):
        components.append("fatty acid")

    # Check for modified amino acids (like DOPA)
    dopa_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H](Cc1ccc(O)c(O)c1)[CX3](=[OX1])[O,N]")
    if mol.HasSubstructMatch(dopa_pattern):
        components.append("modified amino acid (DOPA)")

    # Check for IAN groups
    ian_pattern = Chem.MolFromSmarts("C(#N)Cc1c[nH]c2ccccc12")
    if mol.HasSubstructMatch(ian_pattern):
        components.append("IAN group")

    # Check for CoA
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[CH]([OH])C(C)(C)COP([OH])(=O)OP([OH])(=O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP([OH])([OH])=O")
    if mol.HasSubstructMatch(coa_pattern):
        components.append("Coenzyme A")

    # Check for glutathione core (count as one component if not part of larger peptide)
    gsh_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4H](CCC(=O)[NX3,NX4+][CX4H](CS)[CX3](=O)[NX3,NX4+][CX4H])[CX3](=O)[O-,OH]")
    if mol.HasSubstructMatch(gsh_pattern) and len(components) == 0:
        components.append("glutathione")

    # Check for other significant modifications
    if mol.HasSubstructMatch(Chem.MolFromSmarts("SC(=N)N")): # Thiourea
        components.append("thiourea modification")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)O")): # Sulfonic acid
        components.append("sulfonic acid modification")

    # Final decision
    unique_components = len(set(components))
    if unique_components >= 2:
        return True, "Bioconjugate containing: " + "; ".join(set(components))
    elif unique_components == 1:
        return False, "Only one component found: " + components[0]
    else:
        return False, "No biological components identified"