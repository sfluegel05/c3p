"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# Predefined L-amino acid templates (including modified ones)
L_AMINO_TEMPLATES = {
    'alanine': 'C[C@H](N)C(=O)O',
    'arginine': 'N[C@@H](CCCNC(=N)N)C(=O)O',
    'asparagine': 'N[C@@H](CC(N)=O)C(=O)O',
    'aspartic_acid': 'N[C@@H](CC(=O)O)C(=O)O',
    'cysteine': 'N[C@@H](CS)C(=O)O',
    'glutamic_acid': 'N[C@@H](CCC(=O)O)C(=O)O',
    'glutamine': 'N[C@@H](CCC(N)=O)C(=O)O',
    'glycine': 'NCC(=O)O',
    'histidine': 'N[C@@H](CC1=CN=CN1)C(=O)O',
    'isoleucine': 'CC[C@H](C)[C@H](N)C(=O)O',
    'leucine': 'CC(C)C[C@H](N)C(=O)O',
    'lysine': 'NCCCC[C@H](N)C(=O)O',
    'methionine': 'CSCC[C@H](N)C(=O)O',
    'phenylalanine': 'N[C@@H](CC1=CC=CC=C1)C(=O)O',
    'proline': 'OC(=O)[C@@H]1CCCN1',
    'pyrrolysine': 'C(=O)([C@@H](N)CCCCNC([C@H]1[C@@H](CC=N1)C)=O)O',
    'selenocysteine': 'N[C@@H](C[SeH])C(=O)O',
    'serine': 'N[C@@H](CO)C(=O)O',
    'threonine': 'C[C@H](O)[C@H](N)C(=O)O',
    'tryptophan': 'N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O',
    'tyrosine': 'N[C@@H](CC1=CC=C(O)C=C1)C(=O)O',
    'valine': 'CC(C)[C@H](N)C(=O)O',
    'n-formylmethionine': 'CN(C(=O)SCCC[C@H](N)C(=O)O)C=O',  # Example, adjust as needed
}

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for exactly one amino group and one carboxyl group
    amino_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3&!$(NC=O)]")))
    carboxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1]")))
    
    # Allow N-formyl (pyrrolysine, N-formylmethionine) and other modifications
    if amino_count != 1 and not any('n-formyl' in k for k in L_AMINO_TEMPLATES if 'methionine' in k):
        return False, "Multiple amino groups"
    if carboxyl_count != 1:
        return False, "Multiple carboxyl groups"

    # Check alpha-amino acid pattern with possible N-substituents (like formyl)
    pattern = Chem.MolFromSmarts("[NX3][C@@H]([CX3](=O)[OX2H1])*")
    if not mol.HasSubstructMatch(pattern):
        return False, "No alpha-amino acid backbone with proper stereochemistry"

    # Check against known L-amino acid templates (including modified ones)
    for template_smiles in L_AMINO_TEMPLATES.values():
        template = Chem.MolFromSmiles(template_smiles)
        if template is None:
            continue
        if mol.HasSubstructMatch(template):
            return True, "Matches known proteinogenic amino acid structure"

    # Special case for glycine (non-chiral)
    glycine = Chem.MolFromSmiles(L_AMINO_TEMPLATES['glycine'])
    if mol.HasSubstructMatch(glycine):
        return True, "Glycine structure detected"

    # Check for peptide bonds to exclude dipeptides
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[NX3]")):
        return False, "Contains peptide bond"

    return False, "Does not match any known proteinogenic amino acid structure"