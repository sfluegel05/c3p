"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is classified as a monoamine based on its SMILES string.
    A monoamine contains an amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an aromatic ring (using benzene as a pattern)
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"
    
    # Look for an amino group (-NH2 or primary/secondary/tertiary amine)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")  # this matches primary/secondary/tertiary amines
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amino group found"

    # Look for two-carbon chain connecting amino group to aromatic ring
    chain_pattern = Chem.MolFromSmarts("[$(c)-C-C-N]")  # Aromatic carbon connected to C-C-N
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No C-C chain connecting amino to aromatic ring found"

    return True, "Contains amino group connected to aromatic ring by a two-carbon chain"

# Examples of monoamines, you can test with these SMILES strings
# is_monoamine("CNC[C@H](O)c1ccc(O)c(O)c1")
# is_monoamine("CNCC(O)C1=CC(O)=C(O)C=C1")