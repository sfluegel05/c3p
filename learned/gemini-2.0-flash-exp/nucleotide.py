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
    
    # Define patterns
    # flexible pentose ring
    pentose_pattern1 = Chem.MolFromSmarts("[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]") # 5 C
    pentose_pattern2 = Chem.MolFromSmarts("[C;X4]~[C;X4]~[C;X4]~[O;X2]~[C;X4]") # 4 C and one O

    #  N-glycosidic linkage. The anomeric carbon (C1') must be bonded to an O and a N.
    N_glycosidic_pattern = Chem.MolFromSmarts("[N;!H0]1[#6]~[C;X4]([OX2])~[C;X4]~[C;X4]~[C;X4]~[#6]1")
    N_glycosidic_pattern2 = Chem.MolFromSmarts("[N;!H0]1[#6]~[C;X4]([OX2])~[C;X4]~[C;X4]~[O;X2]~[#6]1")
    
    # Phosphate connected to a sugar oxygen (3' or 5'), and allow for mono, di, and tri phosphates.
    phosphate_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])~[OX2]-[#6]")
    phosphate_pattern2 = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])~[OX2]-[P](=[OX1])([OX2])([OX2])~[OX2]-[#6]")
    phosphate_pattern3 = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])~[OX2]-[P](=[OX1])([OX2])([OX2])~[OX2]-[P](=[OX1])([OX2])([OX2])~[OX2]-[#6]")

    # Check for phosphate group
    if not (mol.HasSubstructMatch(phosphate_pattern) or 
            mol.HasSubstructMatch(phosphate_pattern2) or
            mol.HasSubstructMatch(phosphate_pattern3)
            ):
            return False, "No phosphate group connected to sugar found"
    
    # Check for a pentose sugar
    if not (mol.HasSubstructMatch(pentose_pattern1) or mol.HasSubstructMatch(pentose_pattern2)):
        return False, "No pentose sugar ring found"

    # Check for a nucleobase - an N atom in an aromatic cycle
    nucleobase_pattern = Chem.MolFromSmarts("[n;r]")
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Check for N-glycosidic bond
    if not (mol.HasSubstructMatch(N_glycosidic_pattern) or mol.HasSubstructMatch(N_glycosidic_pattern2)):
       return False, "No N-glycosidic linkage"
    
    
    return True, "Contains a nucleobase, a pentose sugar and a phosphate group connected correctly."