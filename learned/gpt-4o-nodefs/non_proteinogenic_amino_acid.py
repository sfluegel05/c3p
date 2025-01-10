"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # General pattern for amino acids with flexibility for different chain lengths
    amino_acid_pattern = Chem.MolFromSmarts("N[C@]?H([C,R0,N0])C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Lacks recognizable amino acid backbone"

    # Expand uncommon patterns to cover broader set
    uncommon_patterns = {
        "halogenated_side_chain": Chem.MolFromSmarts("[CX4][F,Cl,Br,I]"),
        "nitrile_group": Chem.MolFromSmarts("[CX2]#N"),
        "selenium": Chem.MolFromSmarts("[Se]"),
        "amide_or_imide": Chem.MolFromSmarts("[CX3](=O)[NX3]"),
        "keto_acid": Chem.MolFromSmarts("OC(=O)[CX3](=O)"),
        "proline_like_ring": Chem.MolFromSmarts("C1[NH,C@H](C1)[CX3](=O)O"),
        "beta_hydroxy_acid": Chem.MolFromSmarts("O[C@H](C)[CX3](=O)[OX2H]"),
        "nitro_group": Chem.MolFromSmarts("[NX3](=O)=O"),
        "thioether_linkage": Chem.MolFromSmarts("[SX2]C[SX2]")
    }

    # Check if any of the uncommon patterns match
    for name, pattern in uncommon_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains uncommon feature: {name}"

    # Enhanced stereochemistry check
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    unusual_stereochemistry = False
    for center, stereochemistry in chiral_centers:
        if stereochemistry in ('S', 'R') and mol.GetAtomWithIdx(center).GetAtomicNum() == 6:
            unusual_stereochemistry = True

    if unusual_stereochemistry:
        return True, "Unusual stereochemistry identified (potential D-enantiomer or non-standard configs)"

    return False, "No specific non-proteinogenic features detected"