"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amine group (nitrogen connected to carbon)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary or secondary amine
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Look for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a")  # Any aromatic atom
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic group found"
    
    # Check for alkyl chain linking amine and aromatic groups
    aliphatic_chain_connected = False
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetIsAromatic() != end_atom.GetIsAromatic():
            # Ensure one atom is aromatic and the other is not, suggesting a linking chain
            if (begin_atom.GetAtomicNum() == 7 or end_atom.GetAtomicNum() == 7):  # Verify a nitrogen is involved
                aliphatic_chain_connected = True
                break
    
    if not aliphatic_chain_connected:
        return False, "No alkyl chain linking amine and aromatic group"
    
    return True, "Contains an amine group linked to an aromatic group by an alkyl chain"