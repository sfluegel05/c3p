"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for the basic amino acid backbone structure: an alpha-amino acid chain
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4]([CX4H3,CH2])[CX3](=O)[OX1H1]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Lacks typical amino acid backbone elements"

    # Define expanded uncommon patterns for non-proteinogenic features
    uncommon_patterns = [
        Chem.MolFromSmarts("[CX4]([F,Cl,Br,I])"),         # Halogenated side chains
        Chem.MolFromSmarts("[CX4]#N"),                   # Nitrile group
        Chem.MolFromSmarts("[Se]"),                      # Selenium
        Chem.MolFromSmarts("[CX3](=[OX1])[NX3]"),        # Amide or imide linkages
        Chem.MolFromSmarts("[CX3](=O)[CX3](=O)"),        # Keto-acid group
        Chem.MolFromSmarts("C=O[C@H]1CNC1"),             # Proline-like rings with additional oxygen
        Chem.MolFromSmarts("[OX2CR]([OX1H1])[#6]"),      # Beta-hydroxy acids
        Chem.MolFromSmarts("[CX3]=[N+][O-]"),             # Nitro group
        Chem.MolFromSmarts("[S,C]C[S,C]"),                # Thioether linkages often found in modified methionine
        # Expand as needed to capture other distinctive features
    ]

    for pattern in uncommon_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unusual modification or uncommon side chain indicative of non-proteinogenic amino acids"

    # Check for non-standard amino acid configurations, like D-enantiomers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    d_stereochemistry = False
    for center, stereochemistry in chiral_centers:
        if stereochemistry in ('S', 'R') and mol.GetAtomWithIdx(center).GetAtomicNum() == 6:  # Check for chiral centers in carbon atoms
            d_stereochemistry = True

    if d_stereochemistry:
        return True, "Contains chiral centers with unusual stereochemistry"

    return False, "Does not display uncommon characteristics of non-proteinogenic amino acids"