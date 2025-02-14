"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5
    hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for a pentose sugar
    pentose_pattern = Chem.MolFromSmarts("C1OC(C)C(C)C1")  # basic pentose ring
    if not mol.HasSubstructMatch(pentose_pattern):
       return False, "No pentose sugar ring found"


    # Check for a nucleobase - simplified as N atom connected to a ring system
    nucleobase_pattern = Chem.MolFromSmarts("[N;!H0]~[#6]") # N not connected to H, connected to carbon
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found (simplified)"
    
    # Check if a phosphate is connected to a pentose carbon and another oxygen (sugar)
    phosphate_connection_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])~[#6][OX2]")
    if not mol.HasSubstructMatch(phosphate_connection_pattern):
       return False, "Phosphate not connected to sugar"

    # check for N-glycosidic bond to the pentose
    N_glycosidic_pattern = Chem.MolFromSmarts("[N;!H0]~[#6]1O[#6][#6][#6]1") # base N attached to pentose
    if not mol.HasSubstructMatch(N_glycosidic_pattern):
        return False, "No N-glycosidic linkage"
    
    return True, "Contains a nucleobase, a pentose sugar and a phosphate group connected correctly."