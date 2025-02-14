"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine consists of a sphinganine backbone with a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for sphinganine backbone and N-acyl group (amide bond)
    # the fatty acyl chain can be of any length
    sphinganine_acyl_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CHX4]([CH2X4]O)N[CX3](=[OX1])-[CX4,CX3]")
    if not mol.HasSubstructMatch(sphinganine_acyl_pattern):
        return False, "No N-acylsphinganine substructure found (sphinganine directly connected to acyl group)"

    # Count Carbons - sphinganine has 18 and acyl has at least 2 (acyl C=O and CH2)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 + 2: # minimum carbons for sphinganine + 2 carbons in acyl group
        return False, "Too few carbons for a N-acylsphinganine"
    
    # Count Oxygens - sphinganine has 2, amide has 1, and acyl has at least 1
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3: # minimum oxygens for sphinganine and amide
        return False, "Too few oxygens for a N-acylsphinganine"

    # check molecular weight - N-acylsphinganine is relatively large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
       return False, "Molecular weight too low for N-acylsphinganine"


    return True, "Contains sphinganine backbone with a fatty acyl group attached to the nitrogen"