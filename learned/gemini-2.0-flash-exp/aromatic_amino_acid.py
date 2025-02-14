"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid contains an aromatic ring directly or indirectly connected to the alpha carbon of an amino acid moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring connected via a linker to the alpha carbon of an amino acid
    # [c,n] - aromatic carbon or nitrogen
    # [CX4,CX3]~[CX4,CX3]~* - linker of at least two carbons, followed by an atom
    # [C]([NX3,NX4+,N]) - alpha carbon connected to an amino or imino group
    # C(=O)[O,OH] - carboxylic acid group (neutral or anionic)
    aromatic_amino_acid_linker_pattern = Chem.MolFromSmarts("[c,n]~[CX4,CX3]~[CX4,CX3]~*[C]([NX3,NX4+,N])C(=O)[O,OH]")

    if mol.HasSubstructMatch(aromatic_amino_acid_linker_pattern):
        return True, "Contains an aromatic ring connected via a linker to the alpha carbon of an amino acid moiety"
    
    # Check for aromatic ring directly bonded to alpha carbon of an amino acid moiety
    # This pattern is similar to the linker pattern, but with one single bond between aromatic ring and alpha carbon: [c]~[C]([NX3,NX4+,N])C(=O)[O,OH]
    aromatic_alpha_amino_acid_pattern = Chem.MolFromSmarts("[c,n]~[C]([NX3,NX4+,N])C(=O)[O,OH]")
    if mol.HasSubstructMatch(aromatic_alpha_amino_acid_pattern):
        return True, "Contains an aromatic ring directly bonded to the alpha carbon of an amino acid moiety"

    return False, "No aromatic amino acid moiety found"