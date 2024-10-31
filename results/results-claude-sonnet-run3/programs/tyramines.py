from rdkit import Chem
from rdkit.Chem import AllChem

def is_tyramines(smiles: str):
    """
    Determines if a molecule is a tyramine derivative.
    Tyramines contain a phenol ring connected to an aminoethyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tyramine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for tyramine skeleton:
    # - Aromatic ring with OH group (phenol)
    # - Connected to a chain with an amino group
    tyramine_pattern = Chem.MolFromSmarts('c1([OH])ccc(CC[NX3])cc1')
    meta_tyramine_pattern = Chem.MolFromSmarts('c1c([OH])cc(CC[NX3])cc1')
    
    # Also match variations where the amino group is substituted
    n_acyl_pattern = Chem.MolFromSmarts('c1([OH])ccc(CC[NH]C(=O)*)cc1')
    meta_n_acyl_pattern = Chem.MolFromSmarts('c1c([OH])cc(CC[NH]C(=O)*)cc1')
    
    # Match pattern where there's an OH on the chain
    hydroxy_pattern = Chem.MolFromSmarts('c1([OH])ccc(C[CH]([OH])[NX3])cc1')
    meta_hydroxy_pattern = Chem.MolFromSmarts('c1c([OH])cc(C[CH]([OH])[NX3])cc1')

    if (mol.HasSubstructMatch(tyramine_pattern) or 
        mol.HasSubstructMatch(meta_tyramine_pattern) or
        mol.HasSubstructMatch(n_acyl_pattern) or
        mol.HasSubstructMatch(meta_n_acyl_pattern) or
        mol.HasSubstructMatch(hydroxy_pattern) or
        mol.HasSubstructMatch(meta_hydroxy_pattern)):
        return True, "Contains tyramine skeleton (phenol ring with aminoethyl chain)"

    return False, "Does not contain tyramine skeleton"
# Pr=1.0
# Recall=1.0