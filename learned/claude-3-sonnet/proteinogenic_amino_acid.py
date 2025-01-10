"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get basic formula without isotopes
    mol_no_isotope = Chem.MolFromSmiles(smiles.replace("[2H]", "H").replace("[13C]", "C").replace("[15N]", "N"))
    if mol_no_isotope is None:
        return False, "Invalid SMILES after isotope removal"
        
    # Basic amino acid pattern (non-proline): NH2-CH-COOH
    basic_aa_pattern = Chem.MolFromSmarts("[NX3H2][CX4][CX3](=[OX1])[OX2H1]")
    # Proline pattern: NH-CH-COOH with ring
    proline_pattern = Chem.MolFromSmarts("[NX3H1][CX4][CX3](=[OX1])[OX2H1]")
    
    is_basic_aa = mol_no_isotope.HasSubstructMatch(basic_aa_pattern)
    is_proline = mol_no_isotope.HasSubstructMatch(proline_pattern)
    
    if not (is_basic_aa or is_proline):
        return False, "Missing basic amino acid structure"

    # Check chirality for non-glycine
    formula = CalcMolFormula(mol_no_isotope)
    is_glycine = False
    
    # Special check for glycine and its isotope variants
    if "C2H5NO2" in formula.replace("1", ""):  # Remove isotope numbers
        glycine_pattern = Chem.MolFromSmarts("[NX3H2]CC(=O)O")
        if mol_no_isotope.HasSubstructMatch(glycine_pattern):
            is_glycine = True
    
    if not is_glycine:
        # Check for L configuration in non-glycine amino acids
        chiral_centers = Chem.FindMolChiralCenters(mol_no_isotope, includeUnassigned=True)
        if not chiral_centers:
            return False, "Missing required chirality"
        
        # In SMILES, L-amino acids typically have S configuration at the alpha carbon
        if not any(center[1] == 'S' for center in chiral_centers):
            return False, "Not L configuration"

    # Patterns for specific side chains (simplified)
    side_chains = {
        'alanine': '[CH3]',
        'valine': '[CH](C)C',
        'leucine': '[CH2]C(C)C',
        'isoleucine': '[CH](CC)C',
        'proline': '[CH2][CH2][CH2]N',
        'methionine': '[CH2][CH2]SC',
        'phenylalanine': '[CH2]c1ccccc1',
        'tryptophan': '[CH2]c1c[nH]c2ccccc12',
        'serine': '[CH2]O',
        'threonine': '[CH](O)C',
        'cysteine': '[CH2]S',
        'tyrosine': '[CH2]c1ccc(O)cc1',
        'asparagine': '[CH2]C(N)=O',
        'glutamine': '[CH2][CH2]C(N)=O',
        'aspartic acid': '[CH2]C(O)=O',
        'glutamic acid': '[CH2][CH2]C(O)=O',
        'lysine': '[CH2][CH2][CH2][CH2]N',
        'arginine': '[CH2][CH2][CH2]NC(N)=N',
        'histidine': '[CH2]c1c[nH]cn1',
        'selenocysteine': '[CH2][SeH]',
        'pyrrolysine': '[CH2][CH2][CH2][CH2]NC(=O)[CH]1[CH][CH]N=C1C'
    }

    # Check if molecule matches any of the proteinogenic amino acid patterns
    matched_any = False
    for aa, pattern in side_chains.items():
        aa_pattern = Chem.MolFromSmarts(pattern)
        if aa_pattern and mol_no_isotope.HasSubstructMatch(aa_pattern):
            matched_any = True
            break

    if not (matched_any or is_glycine):
        return False, "Side chain doesn't match any proteinogenic amino acid"

    return True, "Matches proteinogenic amino acid structure with correct configuration"