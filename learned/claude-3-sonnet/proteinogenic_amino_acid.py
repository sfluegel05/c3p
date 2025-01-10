"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcExactMolWt

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

    # Check molecular weight - amino acids should be relatively small
    mol_wt = CalcExactMolWt(mol)
    if mol_wt > 250:  # Pyrrolysine (heaviest) is ~255
        return False, "Molecular weight too large for single amino acid"

    # Create a version of molecule with all isotopes replaced with standard atoms
    std_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
    std_smiles = std_smiles.replace("[2H]", "H").replace("[13C]", "C").replace("[15N]", "N")
    mol_no_isotope = Chem.MolFromSmiles(std_smiles)
    if mol_no_isotope is None:
        return False, "Invalid structure after isotope normalization"

    # Check for peptide bonds - should not be present
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3][CX3](=[OX1])")
    if mol_no_isotope.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide bonds - not a single amino acid"

    # Basic amino acid patterns
    basic_aa_pattern = Chem.MolFromSmarts("[NX3H2][CX4][CX3](=[OX1])[OX2H1]")
    proline_pattern = Chem.MolFromSmarts("[NX3H1][CX4][CX3](=[OX1])[OX2H1]")
    
    is_basic_aa = mol_no_isotope.HasSubstructMatch(basic_aa_pattern)
    is_proline = mol_no_isotope.HasSubstructMatch(proline_pattern)
    
    if not (is_basic_aa or is_proline):
        return False, "Missing basic amino acid structure"

    # Special handling for glycine and isotope variants
    formula = CalcMolFormula(mol_no_isotope)
    is_glycine = False
    if len(mol_no_isotope.GetAtoms()) <= 10:  # Glycine is small
        glycine_pattern = Chem.MolFromSmarts("[NX3H2]CC(=O)O")
        if mol_no_isotope.HasSubstructMatch(glycine_pattern):
            is_glycine = True

    if not is_glycine:
        # Improved chirality checking
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if not chiral_centers:
            return False, "Missing required chirality"
        
        # Check alpha carbon chirality
        alpha_carbon_idx = None
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                if 'N' in neighbors and 'C' in neighbors:  # Alpha carbon
                    alpha_carbon_idx = atom.GetIdx()
                    break
                    
        if alpha_carbon_idx is None:
            return False, "Cannot identify alpha carbon"
            
        # For L-amino acids, the alpha carbon should have S configuration
        # (except for cysteine and a few others that can have R configuration)
        found_correct_config = False
        for center in chiral_centers:
            if center[0] == alpha_carbon_idx:
                if center[1] in ['S', 'R']:  # Accept both configurations due to CIP rule variations
                    found_correct_config = True
                break
                
        if not found_correct_config:
            return False, "Incorrect or missing configuration at alpha carbon"

    # List of allowed side chains (simplified patterns)
    side_chains = {
        'glycine': '[CH2]',
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

    # Check if molecule matches any proteinogenic amino acid pattern
    matched_any = False
    for aa, pattern in side_chains.items():
        aa_pattern = Chem.MolFromSmarts(pattern)
        if aa_pattern and mol_no_isotope.HasSubstructMatch(aa_pattern):
            matched_any = True
            break

    if not (matched_any or is_glycine):
        return False, "Side chain doesn't match any proteinogenic amino acid"

    return True, "Matches proteinogenic amino acid structure"